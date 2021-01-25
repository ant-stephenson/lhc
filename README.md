To build the package:  
Create roxygen documentation  
`devtools::document(roclets = c('rd', 'collate', 'namespace'))`

Convert Rd files into combined pdf (in command line, not R)  
`R CMD Rd2pdf --pdf --title='LHC Package Documentation' -o doc/lhc_package_documentation.pdf man/*.Rd`  

On linux may need to install the following for the above.
`sudo apt-get install texinfo`
`sudo apt-get install texlive-fonts-extra`

Run tests:  
`devtools::test()`  

To load the package:  
`devtools::install_github("ant-stephenson/lhc")`  
`library(lhc)`  


