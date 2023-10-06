# Evaluate the leaf angle probability density function (PDF) from the beta distribution
# using the mean and standard deviation of the leaf inclination angle.
# The leaf angle PDF is calculated for 9 angle classes between 5 and 85 degrees in increments of 10 degrees.

# The variable "leaf" defines the leaf angle distribution type with parameters:
# lad_ave - mean leaf angle (radians)
# lad_std - standard deviation of leaf angle (radians)

leaf <- 'Planophile'
# leaf <- 'Erectophile'
# leaf <- 'Plagiophile'
# leaf <- 'Uniform'
# leaf <- 'Spherical'

switch(leaf,
       'Planophile' = {
         lad_ave <- 26.76 * (pi/180)
         lad_std <- 18.5068 * (pi/180)
       },
       'Erectophile' = {
         lad_ave <- 63.24 * (pi/180)
         lad_std <- 18.4960 * (pi/180)
       },
       'Plagiophile' = {
         lad_ave <- 45.00 * (pi/180)
         lad_std <- 16.2681 * (pi/180)
       },
       'Uniform' = {
         lad_ave <- 45.00 * (pi/180)
         lad_std <- 25.9808 * (pi/180)
       },
       'Spherical' = {
         lad_ave <- 57.30 * (pi/180)
         lad_std <- 21.5485 * (pi/180)
       }
)

# Convert these to the p,q parameters for the beta distribution
num <- 1 - (lad_std*lad_std + lad_ave*lad_ave) / (lad_ave * pi / 2)
den <- (lad_std*lad_std + lad_ave*lad_ave) / (lad_ave*lad_ave) - 1
p <- num / den
q <- ((pi/2) / lad_ave - 1) * p

# Calculate leaf inclination angle probability density function (PDF) and
# fractional abundance (lad) for 9 10-degree bins
dangle <- 10 * (pi/180)                      # Leaf inclination angle increment (radians)
angle <- c(5, 15, 25, 35, 45, 55, 65, 75, 85) # Leaf inclination angle (degrees)
angle <- angle * (pi/180)                    # degrees -> radians

# Initialize vectors to store PDF and lad
beta_pdf <- numeric(length(angle))
beta_lad <- numeric(length(angle))

# Loop through each angle
for (i in 1:length(angle)) {
  x <- angle[i] / (pi/2)
  fp <- x ^ (p - 1)
  fq <- (1 - x) ^ (q - 1)
  beta_pdf[i] <- 2 / pi * fp * fq / beta(p, q)   # Leaf angle PDF 
  beta_lad[i] <- beta_pdf[i] * dangle            # Fraction of leaves in this angle bin
}

# Calculate the known solution
exact_pdf <- numeric(length(angle))
exact_lad <- numeric(length(angle))

# Loop through each angle
for (i in 1:length(angle)) {
  # Exact leaf angle probability density function
  switch(leaf,
         'Planophile' = {
           exact_pdf[i] <- 2 / pi * (1 + cos(2 * angle[i]))
         },
         'Erectophile' = {
           exact_pdf[i] <- 2 / pi * (1 - cos(2 * angle[i]))
         },
         'Plagiophile' = {
           exact_pdf[i] <- 2 / pi * (1 - cos(4 * angle[i]))
         },
         'Uniform' = {
           exact_pdf[i] <- 2 / pi
         },
         'Spherical' = {
           exact_pdf[i] <- sin(angle[i])
         }
  )
  
  # Exact relative leaf angle distribution (fraction)
  exact_lad[i] <- exact_pdf[i] * dangle
}

# Print out fractional abundance and compare with known solution
beta_sum <- sum(beta_lad)
beta_ave <- sum(angle * beta_lad)
exact_sum <- sum(exact_lad)
exact_ave <- sum(angle * exact_lad)

cat("\n")
cat("Leaf type =", leaf, "\n")
cat("     Angle         beta            exact \n")
for (i in 1:length(angle)) {
  cat(sprintf('%10.2f %15.4f %15.4f \n', angle[i]*180/pi, beta_lad[i], exact_lad[i]))
}

cat("\n")
cat("beta distribution \n")
cat("Sum of leaf angle distribution =", beta_sum, "\n")
cat("Mean leaf angle =", beta_ave*180/pi, "\n")
cat("\n")
cat("Exact solution \n")
cat("Sum of leaf angle distribution =", exact_sum, "\n")
cat("Mean leaf angle =", exact_ave*180/pi, "\n")

# Analytical mean leaf angle
fx <- function(x) {
  switch(leaf,
         'Planophile' = {
           x * 2 / pi * (1 + cos(2 * x))
         },
         'Erectophile' = {
           x * 2 / pi * (1 - cos(2 * x))
         },
         'Plagiophile' = {
           x * 2 / pi * (1 - cos(4 * x))
         },
         'Uniform' = {
           x * 2 / pi
         },
         'Spherical' = {
           x * sin(x)
         }
  )
}
analytical_ave <- integrate(fx, 0, pi/2)$value

cat("\n")
cat("Analytical solution \n")
cat("Mean leaf angle =", analytical_ave*180/pi, "\n")

# Calculate Ross index
F1 <- sum(beta_lad[1:3])
F2 <- sum(beta_lad[4:6])
F3 <- sum(beta_lad[7:9])
beta_xl <- 0.5 * (abs(0.134 - F1) + abs(0.366 - F2) + abs(0.5 - F3))
if ((0.5 - F3) < 0) {
  beta_xl <- -beta_xl
}

F1 <- sum(exact_lad[1:3])
F2 <- sum(exact_lad[4:6])
F3 <- sum(exact_lad[7:9])
exact_xl <- 0.5 * (abs(0.134 - F1) + abs(0.366 - F2) + abs(0.5 - F3))
if ((0.5 - F3) < 0) {
  exact_xl <- -exact_xl
}

cat("\n")
cat("Ross index \n")
cat("beta distribution =", beta_xl, "\n")
cat("Exact solution =", exact_xl, "\n")

# Graph PDFs
x <- seq(0, pi/2, length.out = 101)
switch(leaf,
       'Planophile' = {
         y <- 2 / pi * (1 + cos(2 * x))
       },
       'Erectophile' = {
         y <- 2 / pi * (1 - cos(2 * x))
       },
       'Plagiophile' = {
         y <- 2 / pi * (1 - cos(4 * x))
       },
       'Uniform' = {
         y <- 2 / pi
       },
       'Spherical' = {
         y <- sin(x)
       }
)

x <- x * (180/pi)           # radians -> degrees
angle <- angle * (180/pi)   # radians -> degrees

plot(angle, beta_pdf, type = 'l', col = 'blue', xlab = 'Leaf angle (degrees)',
     ylab = 'PDF', main = leaf)
lines(angle, exact_pdf, col = 'red', lty = 2)
lines(x, y, col = 'green', lty = 3)
legend('best', legend = c('beta', 'exact', 'analytical'), col = c('blue', 'red', 'green'), lty = c(1, 2, 3))
