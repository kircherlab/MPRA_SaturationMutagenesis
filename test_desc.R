d <- read.dcf('DESCRIPTION')
cat('DESCRIPTION file parsed successfully!\n')
cat('Package:', d[1,'Package'], '\n')
cat('Version:', d[1,'Version'], '\n')
