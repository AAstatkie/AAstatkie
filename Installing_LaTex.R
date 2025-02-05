#################################
## Installing LaTex in RStudio ##
#################################

# Luckily there is a very nice package that was created for the easy installation
    # of LaTeX in RStudio! Type install.packages("tinytex") into the Console and 
    # press return to run the code.
    # After that is complete, type tinytex::install_tinytex() into the Console 
    # and press return. This should take a bit of time depending on the speed of
    # your computer.
              
    # For some reason, even after a successful installation, sometimes it shows
    # some error/warning messages at the end. Ignore them and check whether it 
    # works here:
                  install.packages("tinytex")
                  library(tinytex)
                  tinytex::install_tinytex()
                  library(tinytex)
                  library(knitr)

    # To check whether it was installed properly
    # 1. Go to File -> New File -> RMarkdown.
    # 2. Then click PDF as the default output option. It will give you example
          # text in the file.
    # 3. Press the Knit button (with the yarn icon) and name the file whatever 
          # you want (Test is always a good option) and put it on your Desktop.
    # 4. It may take a couple of minutes, but you should have a PDF with the same
          # file name (Test.pdf for example) on your Desktop, if it works.
    # 5. If it says Error: LaTeX failed to compile, that means the tinytex
          # installation did not work. Make sure you ran both lines!