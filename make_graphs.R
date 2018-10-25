#!/usr/bin/Rscript --vanilla
suppressMessages(library(tidyverse))

theme_set(theme_minimal())

argv <- commandArgs(trailingOnly=TRUE)

if (length(argv) == 1) {
  n_cells <- Inf
} else if (length(argv) == 2){
  n_cells <- argv[2]
} else {
  stop(call. = FALSE, paste0("Wrong number of arguments.\n\n",
       "  Usage: make_graphs.R <summary-path> [max number of cells]\n\n"))
}

flags <- c("LX", "L1", "L2", "I", "D", "J", "K", "B", "PASS")
short_errors <- factor(flags, rev(flags))

fname1 <- paste0(argv[1], ".errors.csv")
errors <- read_tsv(fname1)
errors <- mutate(errors, Error=factor(Error, rev(flags)))
p1 <- ggplot(errors, aes(Fraction, Error, color=Error, label=round(Fraction, 2))) +
  geom_point() +
  geom_segment(aes(x=0, xend=Fraction, y=Error, yend=Error)) +
  xlim(0, 1) +
  labs(title="Errors in barcode identification") +
  theme(
        axis.title.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size=.25),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")

fname2 <- paste0(argv[1], ".cell_barcodes.csv")
barcodes <- read_tsv(fname2, col_type = "cid") %>%
    mutate(Rank = row_number())

barcodes_df <- bind_rows(barcodes) %>%
  filter(Rank <= n_cells) %>%
  mutate(CumSum = cumsum(Count)) %>%
  mutate(CumFrac = CumSum/max(CumSum))

p2 <- ggplot(barcodes_df, aes(Rank, Count)) +
  geom_line() +
  labs(title="Absolute Read Count", x="Ranked Cells", y=expression("Read Count Per Cell")) +
  scale_x_continuous(
                     labels=function(x) trimws(format(x, big.mark = ","))) +
  scale_y_continuous(trans="log10",
                     labels=function(x) trimws(format(x, big.mark = ","))) +
  theme(
        axis.text.y = element_text(angle=90, hjust=.5),
        panel.grid.major = element_line(size=.25),
        panel.grid.minor = element_blank(),
        legend.position = "none")

p3 <- ggplot(barcodes_df, aes(Rank, CumFrac)) +
  geom_line()+
  labs(title="Cumulative Fraction", x="Ranked Cells", y="Cumulative Fraction of Reads per Cell") +
  scale_x_continuous(
                     labels=function(x) trimws(format(x, big.mark = ","))) +
  theme(
        panel.grid.major = element_line(size=.25),
        panel.grid.minor = element_blank(),
        legend.position = "none")

pdf(paste0(argv[1], ".pdf"), width = 4*1.61803, height = 4)
print(p1)
print(p2)
print(p3)
dev.off()
