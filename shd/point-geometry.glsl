#version 400

// take in a point and output a triangle strip with 4 vertices (aka a "quad")
//
layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

uniform mat4 al_ProjectionMatrix;
uniform float pointSize;

in Vertex {
  vec4 color;
  float size;
}
vertex[];

out Fragment {
  vec4 color;
  vec2 mapping;
}
fragment;

void main() {
  mat4 m = al_ProjectionMatrix;   // rename to make lines shorter
  vec4 v = gl_in[0].gl_Position;  // al_ModelViewMatrix * gl_Position

  float r = pointSize;
  r *= vertex[0].size;

  gl_Position = m * (v + vec4(-r, -r, 0.0, 0.0));
  fragment.color = vertex[0].color;
  fragment.mapping = vec2(-1.0, -1.0);
  EmitVertex();

  gl_Position = m * (v + vec4(r, -r, 0.0, 0.0));
  fragment.color = vertex[0].color;
  fragment.mapping = vec2(1.0, -1.0);
  EmitVertex();

  gl_Position = m * (v + vec4(-r, r, 0.0, 0.0));
  fragment.color = vertex[0].color;
  fragment.mapping = vec2(-1.0, 1.0);
  EmitVertex();

  gl_Position = m * (v + vec4(r, r, 0.0, 0.0));
  fragment.color = vertex[0].color;
  fragment.mapping = vec2(1.0, 1.0);
  EmitVertex();

  EndPrimitive();
}
