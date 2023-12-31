//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================



#include "rtweekend.h"

#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include "moving_sphere.h"
#include "bvh.h"
#include "xy_rect.h"
#include "CornellBox.h"
#include "translate.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#include <omp.h>

#include <iostream>
#include <chrono>
#include <fstream>


color ray_color(const ray& r, const vec3& background, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0)
		return color(0, 0, 0);

    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    // 主要是这里
    // 很好理解，如果depth = 0，也就是在光源处，光最强烈, 并且light的emitted恒定返回光源颜色，color = emitted
    // 而后续的球体和地板的emitted则因为其不会发光，返回的是vec3（0,0,0)，所以需要和不断衰减的光线相加，就这样。
    return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
}


// 读取一下地球贴图！
//hittable_list earth() {
//    int nx, ny, nn;
//    unsigned char* texture_data = stbi_load("input/earthmap.jpg", &nx, &ny, &nn, 0);
//
//    auto earth_surface =
//        make_shared<lambertian>(make_shared<image_texture>(texture_data, nx, ny));
//    auto globe = make_shared<sphere>(vec3(0, 0, 0), 2, earth_surface);
//
//    return hittable_list(globe);
//}

// 矩形光源
hittable_list simple_light() {
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(vec3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(vec3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    auto difflight = make_shared<diffuse_light>(make_shared<constant_texture>(vec3(4, 4, 4)));
    objects.add(make_shared<sphere>(vec3(0, 7, 0), 2, difflight));
    objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));

    return objects;
    //return static_cast<hittable_list>(make_shared<bvh_node>(objects, 0, 1));
}

// 场景
hittable_list random_scene() {
    hittable_list world;

    // 使用纹理
    //auto checker = make_shared<checker_texture>(
    //    make_shared<constant_texture>(vec3(0.2, 0.3, 0.1)),
    //    make_shared<constant_texture>(vec3(0.9, 0.9, 0.9))
    //);
    auto pertext = make_shared<noise_texture>(10);
    auto ground_material = make_shared<lambertian>(pertext);


    world.add(make_shared<sphere>(point3(0,-1000,0), 1000, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    auto texture = make_shared<constant_texture>(albedo);
                    sphere_material = make_shared<lambertian>(texture);
                    // 静止
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));

                    // 运动模糊
					//auto center2 = center + vec3(0, random_double(0, .5), 0);
					//world.add(make_shared<moving_sphere>(
					//	center, center2, 0.0, 1.0, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto texture2 = make_shared<constant_texture>(color(0.4, 0.2, 0.1));
    auto material2 = make_shared<lambertian>(texture2);
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    // 使用BVH！
    return static_cast<hittable_list>(make_shared<bvh_node>(world, 0, 1));
    //return world;
}

hittable_list cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.65, 0.05, 0.05)));
    auto white = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.73, 0.73, 0.73)));
    auto green = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.12, 0.45, 0.15)));
    auto light = make_shared<diffuse_light>(make_shared<constant_texture>(vec3(15, 15, 15)));

    // 外部盒子
    objects.add(make_shared<yz_rect>(0.0, 555, 0.0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0.0, 555, 0.0, 555, 0.0, red));
    objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
    objects.add(make_shared<xz_rect>(0.0, 555, 0.0, 555, 0.0, white));
    objects.add(make_shared<xy_rect>(0.0, 555, 0.0, 555, 555, white));
    objects.add(make_shared<xz_rect>(0.0, 555, 0.0, 555, 555, white));

    // 内部的物体
    shared_ptr<hittable> box1 = make_shared<box>(vec3(0, 0, 0), vec3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    objects.add(box1);

    shared_ptr<hittable> box2 = make_shared<box>(vec3(0, 0, 0), vec3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));
    objects.add(box2);

    return static_cast<hittable_list>(make_shared<bvh_node>(objects, 0, 1));
    //return objects;
}

int main() {

    // Image

    //const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 500;
    //const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int image_height = 500;
    const auto aspect_ratio = double(image_width) / image_height;
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    // World
    //auto world = random_scene();
    auto world = cornell_box();

    // Camera
    vec3 lookfrom(278, 278, -800);
    vec3 lookat(278, 278, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.0;
    auto vfov = 40.0;

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

    // 颜色
    const color background(0, 0, 0);

    // stb需要的图像空间
    unsigned char* data = new unsigned char[image_width * image_height * 3];
	// 颜色指针
	unsigned char* color_ptr = new unsigned char[3];

    // Render
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    // 使用omp多线程框架
    int nthreads, tid;

    // 计时模块
    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel private(tid)
    {
		for (int j = 0; j < image_height; j++) {
			std::cerr << "\rScanlines remaining: " << image_height - j << ' ' << std::flush;
			for (int i = 0; i < image_width; ++i) {
				color pixel_color(0, 0, 0);
				for (int s = 0; s < samples_per_pixel; ++s) {
					auto u = (i + random_double()) / (image_width - 1);
					auto v = (j + random_double()) / (image_height - 1);
					ray r = cam.get_ray(u, v);
					pixel_color += ray_color(r, background, world, max_depth);
				}
                
				write_color(std::cout, color_ptr, pixel_color, samples_per_pixel);
				data[(image_height - j - 1) * image_width * 3 + i * 3] =  color_ptr[0];
				data[(image_height - j - 1) * image_width * 3 + i * 3 + 1] = color_ptr[1];
				data[(image_height - j - 1) * image_width * 3 + i * 3 + 2] = color_ptr[2];
			}
		}
    }
    std::cerr << "\nDone.\n";

    // 结束计时，并将结果写入文件
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::ofstream outFile;
    outFile.open("./output/clock.txt", std::ios::out | std::ios::app);

    std::string place = "公司";
    std::string status = "cornell_box have bvh？";

    outFile << "地点：" << place << std::endl;
    outFile << "情况：" << status << std::endl;
    outFile << "sample: " << samples_per_pixel << std::endl;
    outFile << "本次用时为：" << duration.count() << std::endl;
    outFile << std::endl;
    outFile.close();



    // stride_btye = 一行的比特数
    stbi_write_png("./output/test_bvh.png", image_width, image_height, 3, data, image_width * 3);

    delete[] color_ptr;
    delete[] data;

    return 0;
}
