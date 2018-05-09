#version 330
// for raycasting

in vec3 vertex;
out vec4 vertexOut;

void main()
{
    vertexOut = vec4(vertex, 1.0);
}