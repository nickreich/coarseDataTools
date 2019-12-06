## updates in this edition

 - as requested: fixed improper call to `class(dat) == "data.frame"` with a call to `is(dat, "data.frame")`


## R CMD check results
> On windows-x86_64-devel (r-devel)
  checking sizes of PDF files under 'inst/doc' ... NOTE
  Unable to find GhostScript executable to run checks on size reduction

Unable to figure out whether this is an issue with the server running the checks or something in the package itself that I can/should be changing.

> On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ...NB: need Internet access to use CRAN incoming checks
   NOTE
  Maintainer: ‘Nicholas G. Reich <nick@schoolph.umass.edu>’
  
  Possibly mis-spelled words in DESCRIPTION:
    Lessler (22:5)
    al (20:34, 21:14, 22:16)
    et (20:31, 21:11, 22:13)


These are expected, as they are abbreviations and last names from citations requested by CRAN maintainers.

