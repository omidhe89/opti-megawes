function prepare_for_latex(latexFile, fig, caption)
   fig_name = fig.Name;
    pdf_filename = strcat('.\results\', fig_name, '.pdf');
    latex_dir = strcat('{./results/', fig_name, '.pdf}\n');
    exportgraphics(fig, pdf_filename, 'ContentType', 'image', Resolution=600);
    close(fig);
    fprintf(latexFile, '\\begin{figure}[ht]\n');
    fprintf(latexFile, '\\centering\n');
    fprintf(latexFile, strcat('\\includegraphics[trim={2.9cm 0.15cm 2.9cm 0.15cm}, clip, width=\\textwidth]\n', latex_dir));
    fprintf(latexFile, strcat('\\caption{', caption,'}\n'));
    fprintf(latexFile, strcat('\\label{fig:', fig_name,'}\n'));
    fprintf(latexFile,  '\\end{figure}\n');
end