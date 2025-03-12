# Build a PT form ZERO

1. Sphere::Intersect
    Return root1 > 0 ? root1 : (root2 > 0 ? root2 : 0) to consider the ray starts in the sphere
    Use Vec op = r.o - p; to simplify code

2. camera in main
    Simple implementation. Just emit rays through pixels one by one. No perspective, simple random sampling...
    Notice the format .ppm, which record y first then x.

3. diffuse in radiance