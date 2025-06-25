# Trabalho balanceado
A estrutura dos arquivos está organizada da seguinte forma:
- tba[1..10].ipynb -> código com o solver exato
- tba[1..10]-sa.ipynb -> código com o algoritmo Simulated Annealing puro
- tba[1..10]-hib.ipynb -> código hibrido com o Simulated Annealing servindo de "warm start" do solver puro
- tba[1..10].jl -> código em julia para profiling e visualização do tempo de execução
