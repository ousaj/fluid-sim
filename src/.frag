#version 330 core

in vec3 color;
out vec4 fragColor;
in float fragDrawDisk;

void main() {
	if (fragDrawDisk == 1.0) {
		vec2 p = gl_PointCoord - vec2(0.5);
		if (dot(p, p) > 0.25)
			discard;
	}

	fragColor = vec4(color, 1.0);
}
