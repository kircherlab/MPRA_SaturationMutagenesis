app <- ShinyDriver$new("../", seed = 562346234)
app$snapshotInit("mytest")
app$snapshot()
