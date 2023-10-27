

leaf_area_data = read.csv('../Fortran/02_description_of_ecosystems/leaf_area/leaf_area_density.csv')

plot(leaf_area_data$LAD3, leaf_area_data$Height,  type = 'l', col = 'blue', xlab = 'Leaf area density (m^2 m^{-3})',
     ylab = 'Height (z/h_c)', main = 'Profiles')
lines(leaf_area_data$LAD2, leaf_area_data$Height, col = 'red')
lines(leaf_area_data$LAD1, leaf_area_data$Height, col = 'green')


# % Supplemental program 2.1
# 
# % ----------------------------------------------
#   % Calculate and graph leaf area density profiles
# % ----------------------------------------------
#   

# Parameters for beta distribution
p <- c(2.5, 3.5, 11.5)
q <- c(2.5, 2.0, 3.5)

# Canopy parameters
LAI <- 5
hc <- 10
cat("Leaf area index =", LAI, "\n")

# Create a vector of heights (z) with linearly spaced values
z_min <- 0
z_max <- hc
dz <- 0.1
z <- seq(z_min, z_max, by = dz)

# Initialize lists to store leaf area density profiles
y1 <- numeric(length(z))
y2 <- numeric(length(z))
y3 <- numeric(length(z))

# Calculate the leaf area density profile for each [p,q]
for (i in 1:length(p)) {
  
  sum <- 0
  
  # Loop over each height
  for (j in 1:length(z)) {
    x <- z[j] / hc
    lad <- (LAI / hc) * (x^(p[i]-1) * (1 - x)^(q[i]-1)) / beta(p[i], q[i])
    
    # Numerically sum leaf area for each height
    sum <- sum + lad * dz
    
    # Save output for graphing
    if (i == 1) {
      y1[j] <- lad
    } else if (i == 2) {
      y2[j] <- lad
    } else if (i == 3) {
      y3[j] <- lad
    }
  }
  
  cat("p, q =", p[i], q[i], "\n")
  cat("Leaf area index (numerical) =", sum, "\n")
}

# Make a graph for leaf area density in relation to relative height (z/hc)
z_rel <- z / hc

plot(y3, z_rel, type = 'l', col = 'blue', xlab = 'Leaf area density (m^2 m^{-3})',
     ylab = 'Height (z/h_c)', main = 'Profiles')
lines(y2, z_rel, col = 'red')
lines(y1, z_rel, col = 'green')
