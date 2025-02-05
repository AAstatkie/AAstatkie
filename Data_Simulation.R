# The data simulation practice here is based on the tutorial on https://stirlingcodingclub.github.io/simulating_data/index.html.

# ======================Univariate data simulation==========================
# ======Sampling from a uniform distribution========
# To create a random univariate data of 10 observations with uniform distribution between 0 and 1, use "runif(n, min = 0, max = 1)".

rand_unifs_10<-runif(n=10, min = 0, max = 1) # You will get different numbers each time you run the script.To get the same values each time you run it, set a seed.

# If you want to a large number of observations (e.g., 1000), just increase the value of n.
rand_unifs_1000<-runif(n=1000, min=0, max=1)

# Now, you can summarize the above data (e.g., using histogram).
hist(rand_unifs_1000, col = "magenta", xlab = "Simulated random number", ylab = "Frequency", main = " ", cex.lab=1, cex.axis=1)

# ======Sampling from a normal distribution=========
# To create a random univariate data of 10 observations with uniform/Gaussian distribution with mean 0 and SD 1, use "rnorm()".
rand_norms_10<-rnorm(n=10, mean=0, sd=1)

# To get 1,000 normally distributed values with mean=0 and SD=1:
rand_norms_1000<-rnorm(n=1000, mean=0, sd=1)

# Histogram of the above data:
hist(rand_norms_1000, col = "khaki", xlab = "Random value (X)", main=" ", cex.lab=1, cex.axis=1)

# =====Sampling from a poisson distribution=====
# To create a random univariate data of 10 observations with Poisson distribution with a rate parameter of "lambda", use "rpois()". NOTE: The rate parameter, "lambda", represents both the mean and variance of the distribution, X ~ Poisson(λ).
rand_poissons_10<-rpois(n=10, lambda = 1.5)

# To simulate 1,000 values with poisson distribution with λ=4.5:
rand_poissons_1000<-rpois(n=1000, lambda = 4.5)

# Histogram of the above Poisson-distributed data:
hist(rand_poissons_1000, col="azure", xlab = "Random value (X)", cex.lab=1, cex.axis=1)

# =====Sampling from a binomial distribution=====
# To create a random univariate data of 10 observations with binomial distribution in which each trial generating the binomial outcome is performed only once with a probability of the desired outcome of 0.5, use "rbinom(n=10, size=1, prob=0.5)". 

rand_rbinoms_10<-rbinom(n=10, size = 1, prob = 0.5)

# Too see the output of the above simulation
list(rand_rbinoms_10)

# To simulate 1,000 observations
set.seed(789) # To get the same result each time we run the below code
rand_rbinoms_1000<-rbinom(n=1000, size = 1, prob = 0.5)

# To see the first 10 observations of the above simulation
head(rand_rbinoms_1000, 10)

# To a get table of the frequency distribution of the above simulation
epiDisplay::tab1(rand_rbinoms_1000)

# =====================Random sampling using "sample()"=========================
# Sampling random numbers from a list
# E.g. 1) To select one value at random from a set of values 1-10:
rand_num1<-sample(x = 1:10, size = 1)
rand_num1

# NOTE: In the above function, each vale in the space of 1-10 has equal probability of being selected.

# E.g. 2) If we want to sample 10 values from the values in the range of 1-10, we just set size = 10.
rand_num10<-sample(x = 1:10, size = 10)
rand_num10

# E.g. 3) Sample() by default takes a random sample of values from the range of possible values without replacement. To sample with replacement, we specify the option 'replace = TRUE'.
rand_num10_r<-sample(x = 1:10, size = 10, replace = TRUE)
rand_num10_r

# E.g. 4) We can also specify specific sampling probabilities for each of the values in the range of possible values.For example, to select the values 1-5 with a probability of 0.05 and the values 6-10 with a probability of 0.15, we first create a vector of probabilities, and then specify the vector of probabilities within sample().
prob_vec<-c(rep(x = 0.05, times = 5), rep(x = 0.15, times = 5))
rand_num10_rp<-sample(x = 1:10, size = 10, prob = prob_vec)
rand_num10_rp

# E.g. 5) Sampling random characters from a list
# To create a simulated dataset of n=24 which comprises of three species -- "Species_A, Species_B, and Species_C" -- and we want to do the sampling such that "Species_A" is selected with a probability twice that of "Species_B" and "Species_C":

species <- c("Species_A", "Species_B", "Species_C")
sp_sample <- sample(x = species, size = 24, replace = TRUE, prob = c(0.5, 0.25, 0.25))
print(sp_sample)

# ================Creating a simulated dataframe/dataset====================
# We may create a dataframe of simulated variables.
# E.g.
v1<-runif(n=100, min = 1, max = 5)
v2<-rnorm(n=100, mean = 120, sd = 10)
v3<-rpois(n=100, lambda = 4)
v4<-rbinom(n=100, size = 1, p=0.25)
df<-data.frame(v1, v2, v3, v4)
head(df)

# We may add another simulated variable (v5) to the above dataset (df)
df$v5 <- sample(x = 1:5, size = 100, replace = TRUE)
head(df)

# We can obtain summary statitics of all the variables in df
epiDisplay::summ.data.frame(df)
epiDisplay::tab1(df$v5)

# Histogram of v2 with density curve
hist(df$v2, col = "azure", xlab = "Blood pressure", main = "", probability = TRUE, ylim = c(0, 0.05))
lines(density(df$v2), col="red", lwd=2.0)

# =================***********************************==================

