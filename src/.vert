#version 330 core

layout (location = 0) in vec2 position;
layout (location = 1) in vec3 inColor;

uniform float uPointSize;
uniform vec2 uResolution;
uniform float uDrawDisk;

out vec3 color;
out float fragDrawDisk;

void main() {
	// Converting pixel position to NDC.
	vec2 ndc = (position / uResolution) * 2.0 - 1.0;

	if (uDrawDisk != 1.0) {
		// Offset the center of the point.
		vec2 offset = vec2(uPointSize / 2.0) / uResolution * 2.0;
		offset.y = -offset.y;
		ndc += offset;
	}

	gl_Position = vec4(ndc, 0.0, 1.0);
	color = inColor;
	gl_PointSize = uPointSize;
	fragDrawDisk = uDrawDisk;
}
