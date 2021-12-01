### line2user.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 26 2021 (09:24) 
## Version: 
## Last-Updated: Nov 26 2021 (06:56) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# https://stackoverflow.com/questions/30765866/get-margin-line-locations-in-log-space/30835971#30835971
line2user <- function(line, side, unit) {
    lh <- par('cin')[2] * par('cex') * par('lheight')
    x_off <- diff(grconvertX(c(0, lh), from='inches', to=unit))
    y_off <- diff(grconvertY(c(0, lh), from='inches', to=unit))
    switch(side,
           `1` = grconvertY(-line * y_off, from=unit, to='user'),
           `2` = grconvertX(-line * x_off, from=unit, to='user'),
           `3` = grconvertY(1 + line * y_off, from=unit, to='user'),
           `4` = grconvertX(1 + line * x_off, from=unit, to='user'),
           stop("Argument side must be 1, 2, 3, or 4", call.=FALSE))
}

#----------------------------------------------------------------------
### line2user.R ends here
