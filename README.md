# Build a PT form ZERO

1. Sphere::Intersect
    Return root1 > 0 ? root1 : (root2 > 0 ? root2 : 0) to consider the ray starts in the sphere

    Use Vec op = r.o - p; to simplify code

2. specular in radiance
    Why just using mult for integrate? -- just one direction of the 'input light' so we can ignore integrate!

3. diffuse in radiance
    sample in hemisphere and thus ignore integrate -- no, actually the original smallpt is wrong, we need to divide pdf

4. camera
    learn it from 'ray tracing in one weekend'. The matrix on pbrt do not change the camera ray coordinate! So just sample on position on film plane, then link it with camera. And finally we got a perspective camera. But there are sth wrong in orthographic camera, need to figure out later.

5. refraction
    At first, I suppose the light path of refraction is reverse because we omit light from camera to light source. Actually, we can still regard it as forward light path. We can use classification discussion to prove it: if light goes from air to glass in our camera view, then from light view or 'real view', this light goes from glass to air and will cause total reflection. However, because we know there is definitely an out going light from glass to air('real view'), so we don't need to care about total reflection in this situation one. The second situation is that light goes from glass to air in camera perspective, and is reversed from air to glass in 'real view'. On the contrary, in this situation we need to consider total reflection because we don't sure whether there is still a refractive light going out from glass to air(camera view). All in all, we can regard every light in path tracing as a forward light because of the reversibility.

6. photon mapping, BVH, model input