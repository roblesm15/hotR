
# 1. Create the package structure locally
usethis::create_package("hotR")

# 2. Initialize git inside the package
usethis::use_git()

# 3. Push to GitHub and create the repo simultaneously
usethis::use_github()

# 4. Add a description file

usethis::use_description(fields = list(
  Title = "Structured Threshold Model for Heat-Onset Risk Threshold Estimation",
  Description = "Implements the structured threshold model (STM) for 
    estimating the heat-onset risk threshold (HOT) in mortality time series.",
  `Authors@R` = 'person("Monica", "Robles Fontan", email = "mroblesfontan@g.harvard.edu", 
    role = c("aut", "cre"))'
))

usethis::use_mit_license()

usethis::use_readme_rmd()
