library(RCurl)

user_password <- # "username:password"

url.all88101 <- "ftp://ftp.airtoxics.org/Regional%20Network%20Assessment%20R%20Tools/CorrelationMatrix/all88101.csv"
url.all88502 <- "ftp://ftp.airtoxics.org/Regional%20Network%20Assessment%20R%20Tools/CorrelationMatrix/all88502.csv"
url.allozone <- "ftp://ftp.airtoxics.org/Regional%20Network%20Assessment%20R%20Tools/CorrelationMatrix/ozone_dailysum_2011-2013.csv"
url.ozonemap <- "ftp://ftp.airtoxics.org/Regional%20Network%20Assessment%20R%20Tools/CorrelationMatrix/ozone_map_table.csv"
url.map88101 <- "ftp://ftp.airtoxics.org/Regional%20Network%20Assessment%20R%20Tools/CorrelationMatrix/pm25_map_table_88101.csv"
url.map88502 <- "ftp://ftp.airtoxics.org/Regional%20Network%20Assessment%20R%20Tools/CorrelationMatrix/pm25_map_table_88502.csv"
url.sites <- "ftp://ftp.airtoxics.org/Regional%20Network%20Assessment%20R%20Tools/CorrelationMatrix/small_site_table.csv"
all88101 <- read.csv(text = getURL(url.all88101, userpwd = user_password), stringsAsFactors = FALSE)
all88502 <- read.csv(text = getURL(url.all88502, userpwd = user_password), stringsAsFactors = FALSE)
allozone <- read.csv(text = getURL(url.allozone, userpwd = user_password), stringsAsFactors = FALSE)
ozonemap <- read.csv(text = getURL(url.ozonemap, userpwd = user_password), stringsAsFactors = FALSE)
map88101 <- read.csv(text = getURL(url.map88101, userpwd = user_password), stringsAsFactors = FALSE)
map88502 <- read.csv(text = getURL(url.map88502, userpwd = user_password), stringsAsFactors = FALSE)
sites <- read.csv(text = getURL(url.sites, userpwd = user_password), stringsAsFactors = FALSE)
save(all88101, all88502, allozone, ozonemap, map88101, map88502, sites, file = "removal_bias_data.rdata")
