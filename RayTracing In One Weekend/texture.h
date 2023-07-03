#ifndef TEXTURE_H
#define TEXTURE_H
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

#include "perlin.h"
//#include "stb_image.h"

#include <iostream>

class texture {
public:
    virtual vec3 value(double u, double v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public:
    constant_texture() {}
    constant_texture(vec3 c) : color(c) {}

    virtual vec3 value(double u, double v, const vec3& p) const {
        return color;
    }

public:
    vec3 color;
};

// 棋盘格纹理
class checker_texture : public texture {
public:
    checker_texture() {}
    checker_texture(shared_ptr<texture> t0, shared_ptr<texture> t1) : even(t0), odd(t1) {}

    virtual vec3 value(double u, double v, const vec3& p) const { 
        auto sines = sin(10 * p.x())
                    * sin(10 * p.y())
                    * sin(10 * p.z());

        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }

public:
    shared_ptr<texture> odd;
    shared_ptr<texture> even;
};

// 噪声纹理——大理石
class noise_texture : public texture {
public:
    noise_texture() {}
    noise_texture(double sc) : scale(sc) {}

    virtual vec3 value(double u, double v, const vec3& p) const {
        return vec3(1, 1, 1) * 0.5 * (1.0 + sin(scale * p.z() + 10 * noise.turb(p)));
    }

public:
    perlin noise;
    double scale;
};


//// 图像纹理：读取图像作为纹理！
//class image_texture : public texture 
//{
//public:
//    image_texture() {}
//    image_texture(unsigned char* pixels, int A, int B) 
//        : data(pixels), nx(A), ny(B) {}
//
//    ~image_texture() {
//        delete data;
//    }
//
//    virtual vec3 value(double u, double v, const vec3& p) const {
//        if (data == nullptr)
//            return vec3(0, 1, 1);
//
//        auto i = static_cast<int>((u)*nx);
//        auto j = static_cast<int>((1 - v) * ny - 0.001);
//
//        if (i < 0) i = 0;
//        if (j < 0) j = 0;
//        if (i > nx - 1) i = nx - 1;
//        if (j > ny - 1) j = ny - 1;
//
//        auto r = static_cast<int>(data[3 * i + 3 * nx * j + 0]) / 255.0;
//        auto g = static_cast<int>(data[3 * i + 3 * nx * j + 1]) / 255.0;
//        auto b = static_cast<int>(data[3 * i + 3 * nx * j + 2]) / 255.0;
//
//        return vec3(r, g, b);
//    }
//
//private:
//    unsigned char* data;
//    int nx, ny;
//};
class image_texture : public texture {
public:
    image_texture() {}
    image_texture(unsigned char* pixels, int A, int B)
        : data(pixels), nx(A), ny(B) {}

    ~image_texture() {
        delete data;
    }

    virtual vec3 value(double u, double v, const vec3& p) const {
        // If we have no texture data, then always emit cyan (as a debugging aid).
        if (data == nullptr)
            return vec3(0, 1, 1);

        auto i = static_cast<int>((u)*nx);
        auto j = static_cast<int>((v) * ny - 0.001);

        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i > nx - 1) i = nx - 1;
        if (j > ny - 1) j = ny - 1;

        auto r = static_cast<int>(data[3 * i + 3 * nx * j + 0]) / 255.0;
        auto g = static_cast<int>(data[3 * i + 3 * nx * j + 1]) / 255.0;
        auto b = static_cast<int>(data[3 * i + 3 * nx * j + 2]) / 255.0;

        return vec3(r, g, b);
    }

public:
    unsigned char* data;
    int nx, ny;
};



//class texture  {
//    public:
//        virtual color value(double u, double v, const vec3& p) const = 0;
//};
//
//
//class solid_color : public texture {
//    public:
//        solid_color() {}
//        solid_color(color c) : color_value(c) {}
//
//        solid_color(double red, double green, double blue)
//          : solid_color(color(red,green,blue)) {}
//
//        virtual color value(double u, double v, const vec3& p) const override {
//            return color_value;
//        }
//
//    private:
//        color color_value;
//};
//
//
//class checker_texture : public texture {
//    public:
//        checker_texture() {}
//
//        checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd)
//            : even(_even), odd(_odd) {}
//
//        checker_texture(color c1, color c2)
//            : even(make_shared<solid_color>(c1)) , odd(make_shared<solid_color>(c2)) {}
//
//        virtual color value(double u, double v, const vec3& p) const override {
//            auto sines = sin(10*p.x())*sin(10*p.y())*sin(10*p.z());
//            if (sines < 0)
//                return odd->value(u, v, p);
//            else
//                return even->value(u, v, p);
//        }
//
//    public:
//        shared_ptr<texture> odd;
//        shared_ptr<texture> even;
//};
//
//
//class noise_texture : public texture {
//    public:
//        noise_texture() {}
//        noise_texture(double sc) : scale(sc) {}
//
//        virtual color value(double u, double v, const vec3& p) const override {
//            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
//            // return color(1,1,1)*noise.turb(scale * p);
//            return color(1,1,1)*0.5*(1 + sin(scale*p.z() + 10*noise.turb(p)));
//        }
//
//    public:
//        perlin noise;
//        double scale;
//};
//
//
//class image_texture : public texture {
//    public:
//        const static int bytes_per_pixel = 3;
//
//        image_texture()
//          : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
//
//        image_texture(const char* filename) {
//            auto components_per_pixel = bytes_per_pixel;
//
//            data = stbi_load(
//                filename, &width, &height, &components_per_pixel, components_per_pixel);
//
//            if (!data) {
//                std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
//                width = height = 0;
//            }
//
//            bytes_per_scanline = bytes_per_pixel * width;
//        }
//
//        ~image_texture() {
//            STBI_FREE(data);
//        }
//
//        virtual color value(double u, double v, const vec3& p) const override {
//            // If we have no texture data, then return solid cyan as a debugging aid.
//            if (data == nullptr)
//                return color(0,1,1);
//
//            // Clamp input texture coordinates to [0,1] x [1,0]
//            u = clamp(u, 0.0, 1.0);
//            v = 1.0 - clamp(v, 0.0, 1.0);  // Flip V to image coordinates
//
//            auto i = static_cast<int>(u * width);
//            auto j = static_cast<int>(v * height);
//
//            // Clamp integer mapping, since actual coordinates should be less than 1.0
//            if (i >= width)  i = width-1;
//            if (j >= height) j = height-1;
//
//            const auto color_scale = 1.0 / 255.0;
//            auto pixel = data + j*bytes_per_scanline + i*bytes_per_pixel;
//
//            return color(color_scale*pixel[0], color_scale*pixel[1], color_scale*pixel[2]);
//        }
//
//    private:
//        unsigned char *data;
//        int width, height;
//        int bytes_per_scanline;
//};


#endif
