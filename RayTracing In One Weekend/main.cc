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

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#include <omp.h>

#include <iostream>
#include <chrono>
#include <fstream>


color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    // 这里和阴影无关
	if (depth <= 0)
		return color(0, 0, 0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            //return attenuation * ray_color(scattered, world, depth-1);
            return attenuation;
        return color(0,0,0);
    }

    // 背景
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

// 读取一下地球贴图！
hittable_list earth() {
    int nx, ny, nn;
    unsigned char* texture_data = stbi_load("input/earthmap.jpg", &nx, &ny, &nn, 0);

    auto earth_surface =
        make_shared<lambertian>(make_shared<image_texture>(texture_data, nx, ny));
    auto globe = make_shared<sphere>(vec3(0, 0, 0), 2, earth_surface);

    return hittable_list(globe);
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
}


int main() {

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    // World
    //auto world = random_scene();
    auto world = earth();

    // Camera

    point3 lookfrom(13,2,3);
    point3 lookat(0,0,0);
    vec3 vup(0,1,0);
    //auto dist_to_focus = (lookfrom - lookat).length();
    auto dist_to_focus = 10.0;
    auto aperture = 0.0;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

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
					pixel_color += ray_color(r, world, max_depth);
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
    outFile.open("clock.txt", std::ios::out | std::ios::app);
    outFile << "引入图像贴图" << std::endl;
    outFile << "本次用时为：" << duration.count() << std::endl;
    outFile << "sample: " << samples_per_pixel << std::endl;
    outFile << std::endl;
    outFile.close();



    // stride_btye = 一行的比特数
    stbi_write_png("./output/output.png", image_width, image_height, 3, data, image_width * 3);

    delete[] color_ptr;
    delete[] data;

    return 0;
}
