Hough Transform Explained
===

Simple implementation of the [Hough Transform](https://en.wikipedia.org/wiki/Hough_transform) algorithm that shows with very simple data how the algorithm works in detail, with focus on visualizing what happens. This project helped me a great deal in understanding more about the algorithm and improving also understanding the maths behind it.

Usage (needs [Rust](https://www.rust-lang.org/) installed):

```
$ cd hough-transform-explained
$ cargo run
   Compiling hough-transform-explained v0.1.0
    Finished debug [unoptimized + debuginfo] target(s) in 1.39 secs
     Running `target/debug/hough-transform-explained`
(rho) 127.2628 - round 127
(rho) 127.272415 - round 127
(rho) 129.64806 - round 130
Found line: 254/0 to 0/147 for rho = 127 and theta = 60
(clipped) from 154/58 to 0/147
```

After that you can open `angle-visualization.png`, `houghspace.png` and `detected_lines.png` in the current directory and inspect the output. The article for this code can be found at [my website](http://www.maxbeutel.de/articles/).
