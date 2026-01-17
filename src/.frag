#version 330 core

in vec3 color;
out vec4 fragColor;
in float fragDrawDisk;
in float blurMultiplier;
in float blurParticles;

void main() {
    if (fragDrawDisk == 0.0) { // Grid.
        fragColor = vec4(color, 1.0);
    } else { // Particles.
        vec2 p = gl_PointCoord - vec2(0.5);
        float dist2 = dot(p, p);

        if (dist2 > 0.25)
            discard;

        float alpha = exp(-blurMultiplier * dist2); // Soft edge blur
        if (blurParticles == 1.0) { // Smoothed particles.
            fragColor = vec4(vec3(0.1,0.1,0.1) * alpha, 1.0);
        } else { // Normal particles.
            fragColor = vec4(color * alpha, 1.0);
        }
    }
}
