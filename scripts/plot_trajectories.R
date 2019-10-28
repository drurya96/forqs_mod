#!/usr/bin/env Rscript
#
# plot_trajectories.R
#
# Darren Kessner
# Novembre Lab, UCLA
#


suppressWarnings(require(reshape2, quiet=TRUE))
suppressWarnings(require(ggplot2, quiet=TRUE))


plot_trajectories = function(filename_random, filename_deterministic, title, ylabel, output_filename = "")
{
    if (nchar(output_filename)) pdf(output_filename)

    random = read.table(filename_random)
    deterministic = read.table(filename_deterministic)

    t = cbind(time=1:length(random[,1]), random, deterministic=deterministic$V1) # add time column

    t_melted = melt(t, id.vars="time")
    t_melted = cbind(t_melted, deterministic=(t_melted$variable=="deterministic"))

    a = aes(x=time, y=value, group=variable, color=factor(deterministic))

    g = ggplot(t_melted, a) + 
        geom_line() + 
        scale_color_manual(values=c('grey', 'black'), guide=FALSE) +
        theme_bw() +
        xlab("Generation") + ylab(ylabel) + ggtitle(title)

    print(g)

    if (nchar(output_filename)) dummy=dev.off()
}


main = function()
{
    args = commandArgs(TRUE)
    if (length(args) < 5)
    {
        cat('Usage: plot_trajectories.R filename_random filename_deterministic title ylabel filename_output\n')
        q()
    }

    filename_random = args[1]
    filename_deterministic = args[2]
    title = args[3]
    ylabel = args[4]
    filename_output = args[5]

    plot_trajectories(filename_random, filename_deterministic, title, ylabel, filename_output)
}


main()

