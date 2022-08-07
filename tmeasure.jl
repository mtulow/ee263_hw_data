# This is the script used to produce the data in tomo_data.json.
# Provided for curiosity; YOU DO NOT NEED THIS FILE to solve the problem.

using LinearAlgebra 

npixels = 30    # image size is npixels by npixels
Nd = 35         # number of parallel lines for each angle
Nθ = 35         # number of angles (equally spaced from 0 to pi)
N = Nd * Nθ     # total number of lines (i.e. number of measurements
σ = 0.7         # noise level (standard deviation for normal dist.)

X = # ...the mystery image goes here... (we hid this part!)

x = X[:]    # write the matrix X (the original image) as one big column vector
            # (first column of X goes first, then 2nd column, etc.)

y = zeros(N)
lines_d = zeros(N)
lines_θ = zeros(N)
i = 1

"""
Given an `n`×`n` square grid, and a straight line over that grid `d` units from
the center of the grid (measured in grid lengths orthogonally to the line) at an
angle of `θ` radians from the x-axis, compute length of the line passing through
square in the grid. Return an `n×n` array, whose each element is the length of
the line over the respective grid square.
"""
function line_pixel_length(d, θ, n)

    # for angles in [π/4, 3π/4], flip along diagonal and call recursively
    if π/4 < θ < 3π/4
        L = line_pixel_length(d, π/2-θ, n)
        return collect(transpose(L))
    end

    # for angle in [3π/4, π], redefine line to go in opposite direction
    if θ > π/2
        d = -d
        θ = θ - π
    end

    # for angle in [-π/4, 0], flip along x-axis and call recursively
    if θ < 0
        L = line_pixel_length(-d, -θ, n)
        return reverse(L, dims=1)
    end

    if !(0 ≤ θ ≤ π/2)
        error("invalid angle")
    end

    L = zeros(n, n)
    cosθ = cos(θ)
    sinθ = sin(θ)
    tanθ = sinθ / cosθ

    x0 = n/2 - d * sinθ
    y0 = n/2 + d * cosθ

    y = y0 - x0 .* tanθ
    jy = ceil(Int, y)
    dy = rem(y + n, 1)

    for jx in 1:n
        dynext = dy + tanθ
        if dynext < 1
            if 1 ≤ jy ≤ n
                L[n+1-jy, jx] = 1 / cosθ
            end
            dy = dynext
        else
            if 1 ≤ jy ≤ n
                L[n+1-jy, jx] = (1 - dy) / sinθ
            end
            if 1 ≤ jy + 1 ≤ n
                L[n+1-(jy+1), jx] = (dynext - 1) / sinθ
            end
            dy = dynext - 1
            jy += 1
        end
    end

    return L
end

for iθ = 1:Nθ, id = 1:Nd
    lines_d[i] = 0.7 * npixels * (id - Nd/2 - 0.5) / (Nd/2)
            # equally spaced parallel lines, distance from first to
            # last is about 1.4*n_pixels (to ensure coverage of whole
            # image when theta=pi/4)
    lines_θ[i] = π * iθ / Nθ
            # equally spaced angles from 0 to pi
    L = line_pixel_length(lines_d[i], lines_θ[i], npixels);
            # L is a matrix of the same size as the image
            # with entries giving the length of the line over each pixel
    l = L[:]
            # make matrix L into a vector, as for X
    y[i] = l ⋅ x + σ * randn()
            # l⋅x gives "line integral" of line over image,
            # that is, the intensity of each pixel is multiplied by the
            # length of line over that pixel, and then add for all pixels;
            # a random, Gaussian noise, with std sigma is added to the
            # measurement
    i += 1
            # for next measurement line
end