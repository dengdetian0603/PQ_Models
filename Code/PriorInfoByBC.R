x = rbeta(50000, 1, 1)
y = rbeta(50000, 10, 10)

n.bin = 64
fx = density(x, n = n.bin, from = 0, to = 1)$y + 1e-30
fy = density(y, n = n.bin, from = 0, to = 1)$y + 1e-30

fx = fx/sum(fx) * n.bin
fy = fy/sum(fy) * n.bin
dr = 1/n.bin

BC = sum(sqrt(fx * fy)) * dr
BC

BD = -log2(BC)
BD

KL.xtoy= sum(fx * log2(fx/fy)) * dr
KL.xtoy

Hxy = -sum(fx * log2(fy/n.bin)) * dr
Hxy

Hx = -sum(fx * log2(fx/n.bin)) * dr
Hx

Hx / Hxy

sum(fx * log2(fy)) * dr