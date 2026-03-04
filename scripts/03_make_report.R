cat > scripts/03_make_report.R <<'EOF'
#!/usr/bin/env Rscript
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Missing package: rmarkdown")

rmarkdown::render(
  input = "reports/report.Rmd",
  output_file = "report.html",
  output_dir = "results",
  quiet = TRUE
)
message("Wrote: results/report.html")
EOF

chmod +x scripts/03_make_report.R
