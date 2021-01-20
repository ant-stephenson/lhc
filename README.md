To build the package:  
Create roxygen documentation  
`devtools::document(roclets = c('rd', 'collate', 'namespace'))`

Convert Rd files into combined pdf (in command line, not R):
`R CMD Rd2dvi --pdf --title='Test of foo' -o /tmp/foo.pdf man/*.Rd`

Run tests:  
`devtools::test()`  

To load the package:  
`devtools::install_github("ant-stephenson/lhc")`  
`library(lhc)`  


