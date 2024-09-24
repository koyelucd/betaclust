## betaclust Package Information
This is the first submission of this package to CRAN. A few notes about examples and tests: we use `\donttest{}` for examples because all functions use heavy computation. The example data size has been reduced and the runtime is below 5 seconds. But on checking the package build in ubuntu-gcc-release and fedora-clang-devel the elapsed time was observed to be greater than 5 seconds. Hence the examples have been wrapped using `\donttest{}`.

## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
0 errors, 0 warnings, 1 note
  NOTE: New submission

## Resubmission 2022-12-02
This is the second submission of betaclust to CRAN. Comments from reviewer Victoria Wimmer:
```
If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
Missing Rd-tags:
     plot.betaclust.Rd: \value

All your examples are wrapped in \donttest{} and therefore do not get tested.
Please unwrap the examples if that is feasible and if they can be executed in < 5 sec for each Rd file or create additionally small toy examples to allow automatic testing.

You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object.
Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console.
(except for print, summary, interactive functions)
In R/betaclust.R  cat("fitting ...\n")

Please ensure that you do not use more than 2 cores in your examples, vignettes, etc.

Please fix and resubmit.

Best,
Victoria Wimmer
```
In this submission I have:
  - Written references in the description of DESCRIPTION file in the form: "authors (year) <arXiv:...>".  
  
  - Added \value to plot.betaclust.Rd file.  
  
  - I have replaced cat("fitting ...\n") to message("fitting ...\n") as suggested in betaclust.R file.  
  
  - I have ensured that not more than 2 cores are being used in the examples, vignettes, etc.  
  
  - I have removed `\donttest{}` from some examples which are taking <5 secs.  
  
 
 ## Test environments  
 
- R-hub windows-x86_64-devel (r-devel)  

- R-hub macOS 10.13.6 High Sierra, R-release, brew  

- R-hub Debian Linux, R-release, GCC   

- R-hub ubuntu-gcc-release (r-release)  

- R-hub fedora-clang-devel (r-devel)   


## R CMD check results
### Local check
`devtools::check()` result:  

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` Windows result:  

  Status: OK
  
  0 errors ✔ | 0 warnings ✔ | 0 notes ✔  
  
`rhub::check_on_debian()` Debian Linux result:  

  Status: OK  
  
  0 errors ✔ | 0 warnings ✔ | 0 notes ✔  
  
`rhub::check_on_macos()` macOS result:  
  Status: OK  
  
  0 errors ✔ | 0 warnings ✔ | 0 notes ✔  
  
`rhub::check_rhub()` Ubuntu Linux result:  
* checking CRAN incoming feasibility ... NOTE  
Maintainer: ‘Koyel Majumdar <koyel.majumdar@ucdconnect.ie>’  

New submission  

* checking examples ... [7s/83s] NOTE  
Examples with CPU (user + system) or elapsed time > 5s  
  0 errors ✔ | 0 warnings ✔ | 2 notes ✖  
  
`rhub::check_rhub()` Fedora Linux result:  
* checking CRAN incoming feasibility ... NOTE  
Maintainer: ‘Koyel Majumdar <koyel.majumdar@ucdconnect.ie>’  

New submission  
* checking examples ... [8s/84s] NOTE  
Examples with CPU (user + system) or elapsed time > 5s  
  0 errors ✔ | 0 warnings ✔ | 2 notes ✖  
  
  
Notes:  

  - The first note is due to New submission and can be ignored.  
  
  - The second note in rhub Linux only occurs in the 2 default rhub Linux environments (Ubuntu and Fedora). It doesn't happen on Windows, Debian Linux or Mac builders, and this behavior is not being replicated locally, so this note can be related to rhub and not the code itself.  
  
## Resubmission 2023-10-02
This is the third submission of betaclust to CRAN.
In this submission I have:
  - Added function for AUC and WD metrics calculation for identifying DMCs as per manuscript resubmission corrections.  
 
 ## Test environments  
 
- R-hub windows-x86_64-devel (r-devel)  

- devtools macOS 13.3.1 

- R-hub Debian Linux, R-release, GCC   

- R-hub ubuntu-gcc-release (r-release)  

- R-hub fedora-clang-devel (r-devel)   


## R CMD check results
### Local check
`devtools::check()` result:  

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` Windows result:  

  Status: OK
  
  0 errors ✔ | 0 warnings ✔ | 3 notes ✖
  
  ❯ checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
    
`rhub::check_on_debian()` Debian Linux result:  

  Status: OK  
  
  0 errors ✔ | 0 warnings ✔ | 0 notes ✔  
  
`devtools::check_mac_release()` macOS result:  
  Status: OK  
  
  0 errors ✔ | 0 warnings ✔ | 0 notes ✔  
  
`rhub::check_rhub()` Ubuntu Linux result:  

* checking examples ... [4s/26s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
         user system elapsed
beta_kr 1.238  0.018   9.350
beta_k  0.274  0.008   5.204
beta_kn 0.180  0.004   5.071

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable
  
`rhub::check_rhub()` Fedora Linux result:  

* checking examples ... [4s/24s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
         user system elapsed
beta_kr 1.271  0.018   9.032

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable
  
  
Notes:  
  
  - The notes in rhub only occurs in the rhub Windows environment. They are mentioned in the open issues of the rhub package.
  - The note for both Ubuntu and Fedora ('tidy' and 'V8'), which, as far as I understood, is caused by some problem on the side of the testing platform.
  
  
## Resubmission 2024-09-24
This is the fourth submission of betaclust to CRAN.
In this submission I have:
  - Corrected a code issue while calculating ICL criteria for optimal model selection.  
 
 ## Test environments  
 
- R-hub windows-x86_64-devel (r-devel)  

- devtools macOS 13.3.1 

- R-hub Debian Linux, R-release, GCC   

- R-hub ubuntu-gcc-release (r-release)  

- R-hub fedora-clang-devel (r-devel) 

## R CMD check results
### Local check
`devtools::check()` result:  

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
