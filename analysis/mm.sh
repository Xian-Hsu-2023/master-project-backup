#!/bin/bash
# make short animations for the plot produced by mpi_pythonplot.sh
# please modify the path for your own purpose
export TEST="magtest057"

mkdir animation/${TEST}
export DIR="figure/${TEST}_*/"
alias mf="mv *.mp4 ../../animation/${TEST}/"
MM()
{
    ffmpeg -r 10 -pattern_type glob -i "*.png" -vcodec libx264 -s 2048x2048 -pix_fmt yuv420p $1.mp4
}

alias cb='cd ../'
alias MMM='MM ${PWD##*/}'
for dir in $DIR;
    do
        cd $dir; MMM; mf; cb; cb;
        echo complete animation from $dir
    done

echo Run complete at `date`
