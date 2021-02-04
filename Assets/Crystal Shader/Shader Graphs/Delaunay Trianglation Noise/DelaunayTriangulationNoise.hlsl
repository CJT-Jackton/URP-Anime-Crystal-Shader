#ifndef DELAUNAY_TRIANGULATION_NOISE_INCLUDED
#define DELAUNAY_TRIANGULATION_NOISE_INCLUDED

float2 delaunay_noise_randomVector(float2 UV, float offset)
{
    float2x2 m = float2x2(15.27, 47.63, 99.41, 89.98);
    UV = frac(sin(mul(UV, m)) * 46839.32);
    return float2(sin(UV.y * offset), cos(UV.x * offset)) * 0.25;
}

// Signed distance to a line crossing (p0, p1) segment.
float distLine(float2 p0, float2 p1)
{
    float2 e0 = p1 - p0;
    return dot(p0, normalize(float2(e0.y, -e0.x)));
}

// Use "Parabolic lifting" method to calculate if two triangles are about to flip.
// This is actually more reliable than circumscribed circle method.
// The technique is based on duality between Delaunay Triangulation
// and Convex Hull, where DT is just a boundary of convex hull
// of projected seeds onto paraboloid.
// We project (h1 h2 h3) triangle ontot paraboloid
// and return the distance of the origin
// to a plane crossing projected triangle.
float flipDistance(float2 h1, float2 h2, float2 h3)
{
    // Project triangle on paraboloid.
    float3 g1 = float3(h1.x, h1.y, dot(h1, h1));
    float3 g2 = float3(h2.x, h2.y, dot(h2, h2));
    float3 g3 = float3(h3.x, h3.y, dot(h3, h3));

    return dot(g1, cross(g3 - g1, g2 - g1));
}

// Find the distance to the closest Delaunay edge in quad (h0, h1, h2, h3).
void delaunayQuad(float2 h0, float2 h1, float2 h2, float2 h3, float2 pos0, float2 pos1, float2 pos2, float2 pos3, inout float distance, inout float2 center)
{
    // Get the distance to quad edges.
    // note: in general it can be concave, but we don't care.
    float minDistance = min(min(distLine(h0, h1), distLine(h1, h2)), min(distLine(h2, h3), distLine(h3, h0)));

    if (minDistance < 0.0) return; // Outside the quad

    // Calculate flip distance relative to h2.
    float dc = flipDistance(h0 - h2, h1 - h2, h3 - h2);

    // Flipping rotates diagonal from (h0 h2) to (h3 h1).
    float2 diagonalPosition0 = (dc > 0.0) ? h3 : h0;
    float2 diagonalPosition1 = (dc > 0.0) ? h1 : h2;

    // Calculate distance to diagonal.
    float distToDiagonal = distLine(diagonalPosition0, diagonalPosition1);

    // Calculate distance to diagonal (positive from both sides).
    minDistance = min(minDistance, abs(distToDiagonal));

    distance = min(minDistance, distance);

    if (dc > 0.0)
    {
        if (distToDiagonal > 0.0)
            center = (pos1 + pos2 + pos3) / 3.0;
        else
            center = (pos0 + pos1 + pos3) / 3.0;
    }
    else
    {
        if (distToDiagonal > 0.0)
            center = (pos0 + pos2 + pos3) / 3.0;
        else
            center = (pos0 + pos1 + pos2) / 3.0;
    }
}

// Final function visits 4 quads around the center cell.
void delaunayTriangulation(float2 UV, float AngleOffset, out float minDistance, out float2 center)
{
    float2 grid = floor(UV);

    // Cache sites coordinate.
    float2 site[9];
    float2 h[9];

    // Iterate sites.
    for (int y = -1; y <= 1; y++)
    {
        for (int x = -1; x <= 1; x++)
        {
            float2 lattice = float2(x, y);
            float2 offset = float2(0.5, 0.5) + delaunay_noise_randomVector(grid + lattice, AngleOffset);
            site[(x + 1) + (y + 1) * 3] = grid + lattice + offset;
            h[(x + 1) + (y + 1) * 3] = site[(x + 1) + (y + 1) * 3] - UV;
        }
    }

    minDistance = 8.0;
    center = float2(0.0, 0.0);

    delaunayQuad(h[3], h[0], h[1], h[4], site[3], site[0], site[1], site[4], minDistance, center);
    delaunayQuad(h[4], h[1], h[2], h[5], site[4], site[1], site[2], site[5], minDistance, center);
    delaunayQuad(h[7], h[4], h[5], h[8], site[7], site[4], site[5], site[8], minDistance, center);
    delaunayQuad(h[6], h[3], h[4], h[7], site[6], site[3], site[4], site[7], minDistance, center);
}

// Delaunay Triangluation
// https://www.shadertoy.com/view/4sKyRD
void Unity_Delaunay_Triangulation_float(float2 UV, float AngleOffset, float CellDensity, out float Out, out float2 Cells)
{
    delaunayTriangulation(UV * CellDensity, AngleOffset, Out, Cells);
}

#endif // DELAUNAY_TRIANGULATION_NOISE_INCLUDED