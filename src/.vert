#version 330 core

layout (location = 0) in vec2 position;
layout (location = 1) in vec3 inColor;

uniform float uPointSize;
uniform vec2 uResolution;
uniform float uDrawDisk;
uniform bool uBlurParticles;
uniform float uGlowMultiplier;
uniform float uBlurMultiplier;

out vec3 color;
out float fragDrawDisk;
out float blurMultiplier;
out float blurParticles;

void main() {
	vec2 ndc = (position / uResolution) * 2.0 - 1.0;

	if (uDrawDisk == 0.0) {
		// Offset the center of the point.
		vec2 offset = vec2(uPointSize / 2.0) / uResolution * 2.0;
		offset.y = -offset.y;
		ndc += offset;
	}

	gl_Position = vec4(ndc, 0.0, 1.0);
	color = inColor;
	if (uDrawDisk == 1.0 && uBlurParticles) {
		gl_PointSize = uPointSize * uGlowMultiplier;
	} else {
		gl_PointSize = uPointSize;
	}

	fragDrawDisk = uDrawDisk;
	blurMultiplier = uBlurMultiplier;
	blurParticles = uBlurParticles ? 1.0 : 0.0;
}
