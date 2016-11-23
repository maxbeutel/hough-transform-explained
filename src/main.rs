extern crate nalgebra as na;
extern crate image;

use std::f32;
use std::fs::File;
use std::path::Path;

use image::{ImageBuffer, GenericImage};

// -- actual implementation --
fn hough_transform() {
    let mut img = image::open(&Path::new("initial-image.png")).unwrap();

    let (img_width, img_height) = img.dimensions();

    // calculate the longest line that can occur in the image (hypotenuse)
    let max_line_length = (img_width as f32).hypot(img_height as f32).ceil();

    let mut angle_visualization_img = ImageBuffer::new(img_width, img_height);

    // Create the accumulator matrix, x-axis is angle theta, y-axis is distance rho
    let mut accumulator: na::DMatrix<u32> = na::DMatrix::new_zeros(180, max_line_length as usize);

    // Flip y axis when reading image, as the image has 0/0 at the top left,
    // but it's easier for us to work with when the origin 0/0 is at the bottom left
    for y in 0..img_height {
        for x in 0..img_width {
            let y_flipped = img_height - 1 - y;

            if !detect_edge(x, y_flipped) {
                continue;
            }

            // we found an edge, now try to draw lines going through this point, using
            // the current point as the origin 0,0

            for i in 0..3 {
                let theta = (i as f32 + 1.0) * 30.0;
                //println!("at {}/{} angle {}", x, y_flipped, theta);

                // -- start visualization of the angle, NOT hough transform related --
                // calculate remaining angle in the triangle (alpha is "angle", beta is 90)
                let angle_beta = 90f32;
                //let angle_gamma = 180f32 - (theta + angle_beta);

                // one side of the triangle is know, it's a, opposite of angle alpha, so the "height" of the triangle
                let length_a = x as f32;

                // now use sinus rule to find the remaining sides
                // don't forget to convert to radians before using trigonometry functions
                let length_b = length_a * angle_beta.to_radians().sin() / theta.to_radians().sin();
                //let length_c = length_b * angle_gamma.to_radians().sin() / angle_beta.to_radians().sin();

                // println!("length_b {}", length_b);
                // println!("length_c {}", length_c);

                // now calculate end of line, with length b of the triangle we just used
                // p1_x/p2_y is the origin so at 0/0
                let p2_x_relative = 0f32 + length_b * theta.to_radians().cos();
                let p2_y_relative = 0f32 + length_b * theta.to_radians().sin();

                //println!("(relative) from {}/{} to {}/{}", x, y_flipped, p2_x_relative, p2_y_relative);

                let p2_x = (p2_y_relative - x as f32).round();
                let p2_y = (p2_x_relative + y_flipped as f32).round();

                //println!("(absolute) from {}/{} to {}/{}", x, y_flipped, p2_x, p2_y);

                // clip the line to fit our image canvas
                let mut clipped_x1 = 0.0;
                let mut clipped_y1 = 0.0;
                let mut clipped_x2 = 0.0;
                let mut clipped_y2 = 0.0;

                liang_barsky(
                    0.0, img_width as f32 - 1.0, 0.0, img_height as f32 - 1.0,
                    x as f32, y_flipped as f32, p2_x as f32, p2_y as f32,
                    &mut clipped_x1, &mut clipped_y1, &mut clipped_x2, &mut clipped_y2
                );

                //println!("(clipped) from {}/{} to {}/{}", clipped_x1.round(), clipped_y1.round(), clipped_x2.round(), clipped_y2.round());

                draw_line(
                    &mut angle_visualization_img,
                    clipped_x1.round() as i32,
                    img_height as i32 - clipped_y1.round() as i32,
                    clipped_x2.round() as i32,
                    img_height as i32 - clipped_y2.round() as i32,
                    image::Rgb([255, 255, 255])
                    //(255, 255, 255)
                );
                // -- end visualization of the angle --

                // this is the HOUGH TRANSFROM:
                // calculate rho
                let rho = (x as f32) * theta.to_radians().cos() + (y_flipped as f32) * theta.to_radians().sin();
                let rho_rounded = rho.round();

                // theta at 60 degree is the angle that has all three data points on the same line
                // they will output similar rho values
                if theta == 60.0 { println!("(rho) {} - round {}", rho, rho_rounded); }

                accumulator[(theta as usize, rho_rounded as usize)] += 1;
            }
        }
    }

    let ref mut angle_visualization_fout = File::create(&Path::new("angle-visualization.png")).unwrap();
    let _ = image::ImageRgb8(angle_visualization_img).save(angle_visualization_fout, image::PNG);


    // dump hough space as image
    let houghspace_img_width = accumulator.nrows() as u32;
    let houghspace_img_height = accumulator.ncols() as u32;

    let mut houghspace_img = ImageBuffer::new(houghspace_img_width, houghspace_img_height);

    // this is still very ugly, cloning the vector to get the maximum votes value
    let accu_clone = accumulator.clone().into_vector();
    let max_accumulator_value = *accu_clone.iter().max().unwrap();

    for theta in 0..houghspace_img_width {
        for rho in 0..houghspace_img_height {
            let n = na::min(((accumulator[(theta as usize, rho as usize)] as f32) * 255.0 / (max_accumulator_value as f32)).round() as u32, 255) as u8;
            let pixel = image::Rgb([n, n, n]);

            houghspace_img[(theta, houghspace_img_height - rho - 1)] = pixel;
        }
    }

    let ref mut houghspace_fout = File::create(&Path::new("houghspace.png")).unwrap();
    let _ = image::ImageRgb8(houghspace_img).save(houghspace_fout, image::PNG);


    // find lines in image based on accumulator maximum
    let accumulator_threshold = 2;

    for theta in 0..houghspace_img_width {
        for rho in 0..houghspace_img_height {
            let accumulator_value = accumulator[(theta as usize, rho as usize)];

            if accumulator_value < accumulator_threshold {
                continue;
            }

            let p1_x = rho as f32 / (theta as f32).to_radians().cos();
            let p1_y = 0.0f32;

            let p2_x = 0.0f32;
            let p2_y = rho as f32 / (theta as f32).to_radians().sin();

            println!("Found line: {}/{} to {}/{} for rho = {} and theta = {}", p1_x.round(), p1_y.round(), p2_x.round(), p2_y.round(), rho, theta);

            // clip the line to fit our image canvas,
            // as the resulting p1/p2 might be outside the image dimensions
            let mut clipped_x1 = 0.0;
            let mut clipped_y1 = 0.0;
            let mut clipped_x2 = 0.0;
            let mut clipped_y2 = 0.0;

            liang_barsky(
                0.0, img_width as f32 - 1.0, 0.0, img_height as f32 - 1.0,
                p1_x, p1_y, p2_x, p2_y,
                &mut clipped_x1, &mut clipped_y1, &mut clipped_x2, &mut clipped_y2
            );

            println!("(clipped) from {}/{} to {}/{}", clipped_x1.round(), clipped_y1.round(), clipped_x2.round(), clipped_y2.round());

            draw_line(
                &mut img, // draw line in original image, showing the ine we found
                clipped_x1.round() as i32,
                img_height as i32 - clipped_y1.round() as i32,
                clipped_x2.round() as i32,
                img_height as i32 - clipped_y2.round() as i32,
                image::Rgba([142, 250, 0, 1])
            );
        }
    }

    let ref mut detected_lines_fout = File::create(&Path::new("detected_lines.png")).unwrap();
    let _ = img.save(detected_lines_fout, image::PNG);
}

// -- utility functions --

// Based on http://stackoverflow.com/questions/34440429/draw-a-line-in-a-bitmap-possibly-with-piston
fn draw_line<T: GenericImage>(img: &mut T, x0: i32, y0: i32, x1: i32, y1: i32, pixel: T::Pixel) {

    // Create local variables for moving start point
    let mut x0 = x0;
    let mut y0 = y0;

    // Get absolute x/y offset
    let dx = if x0 > x1 { x0 - x1 } else { x1 - x0 };
    let dy = if y0 > y1 { y0 - y1 } else { y1 - y0 };

    // Get slopes
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };

    // Initialize error
    let mut err = if dx > dy { dx } else {-dy} / 2;
    let mut err2;

    loop {
        // Set pixel
        img.put_pixel(x0 as u32, y0 as u32, pixel);

        // Check end condition
        if x0 == x1 && y0 == y1 { break };

        // Store old error
        err2 = 2 * err;

        // Adjust error and start position
        if err2 > -dx { err -= dy; x0 += sx; }
        if err2 < dy { err += dx; y0 += sy; }
    }
}

// Liang-Barsky function by Daniel White @ http://www.skytopia.com/project/articles/compsci/clipping.html
// I ported this function almost 1:1 from C++, in all it's ugliness
// This is not related to Hough Transform and is only used for visualization of the found lines in the original image
#[allow(unused_assignments)]
fn liang_barsky (
    edge_left: f32, edge_right: f32, edge_bottom: f32, edge_top: f32,   // Define the x/y clipping values for the border.
    x0src: f32, y0src: f32, x1src: f32, y1src: f32,                 // Define the start and end points of the line.
    x0clip: &mut f32, y0clip: &mut f32, x1clip: &mut f32, y1clip: &mut f32 // The output values, so declare these outside.
) -> bool {

    let mut t0: f32 = 0.0; let mut t1: f32 = 1.0;
    let xdelta = x1src - x0src;
    let ydelta = y1src - y0src;
    let mut p = 0.0f32;
    let mut q = 0.0f32;
    let mut r = 0.0f32;

    for edge in 0..4 {   // Traverse through left, right, bottom, top edges.
        if edge == 0 {  p = -xdelta;    q = -(edge_left - x0src);   }
        if edge == 1 {  p = xdelta;     q =  edge_right - x0src;    }
        if edge == 2 {  p = -ydelta;    q = -(edge_bottom - y0src); }
        if edge == 3 {  p = ydelta;     q =  edge_top - y0src;      }
        r = q/p;

        if p == 0.0 && q < 0.0 {
            // Don't draw line at all. (parallel line outside)
            return false;
        }

        if p < 0.0 {
            if r as f32 > t1 {
                // Don't draw line at all.
                return false;
            } else if r as f32 > t0 {
                // Line is clipped!
                t0 = r as f32;
            }
        } else if p > 0.0 {
            if (r as f32) < t0 {
                // Don't draw line at all.
                return false;
            }
            else if (r as f32) < t1 {
                // Line is clipped!
                t1 = r as f32;
            }
        }
    }

    *x0clip = x0src as f32 + t0 as f32 * xdelta as f32;
    *y0clip = y0src as f32 + t0 as f32 * ydelta as f32;
    *x1clip = x0src as f32 + t1 as f32 * xdelta as f32;
    *y1clip = y0src as f32 + t1 as f32 * ydelta as f32;

    return true;
}

// fake edge detection, refer to the article to understand why/how/what
fn detect_edge(x: u32, y: u32) -> bool {
    if x == 64 && y == 110 { return true; }
    if x == 90 && y == 95 { return true; }
    if x == 119 && y == 81 { return true; }
    return false;
}

// -- main --
fn main() {
    hough_transform();
}
