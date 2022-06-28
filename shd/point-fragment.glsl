#version 400

in Fragment {
  vec4 color;
  vec2 mapping;
}
fragment;

layout(location = 0) out vec4 fragmentColor;

void main() {
  float r = dot(fragment.mapping, fragment.mapping);
  if (r > 1) discard;
  fragmentColor = vec4(fragment.color.rgb, 1 - r * r);
}