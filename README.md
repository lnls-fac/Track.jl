# Track.jl
> Informacoes simples ~ por enquanto!
## Instalacao
Para instalar e testar o pacote basta possuir julia na maquina. 

Julia tambem pode ser instalado dentro de um ambiente CONDA/MAMBA (como o mamba-env ```sirius```). 

A instalacao do julia pode ser feita por:
* No linux: ```sudo snap install julia```
* No Conda: ```conda install -c conda-forge julia```
* No Mamba: ```mamba install -c conda-forge julia```
>Para conda/mamba: recomendo a utilizacao da flag ```--freeze-installed``` para nao haver interferencia com os pacotes do ambiente

A instalacao do pacote Track.jl pode ser feita da seguinte maneira:
- Clone o repositório;
- Entre na pasta do repositório (onde estão os arquivos: Track.jl e Project.toml);
- Ative  ```REPL julia``` com o comando ```julia``` no terminal;
- Entre no modo "```pkg```" (gerenciador de pacotes) teclando: "```]```";
- No modo "```pkg```", execute o comando: "```dev .```";
- Saia do modo "```pkg```", com a tecla "Backslash" ($\leftarrow$);
- Saia do  ```REPL julia``` com o comando ```exit()```;
- Ative novamente o ```REPL julia```;
- Execute o comando ```using Track```.
> Espero que tenha funcionado :D

## Funcionalidades

O pacote Track.jl disponibiliza funções como, por exemplo:
* ``` StorageRing.create_accelerator() ``` constrói o modelo do anel de armazenamento do SIRIUS; 
* ``` Pos() ``` Cria uma partícula (vetor 6-D: $[x, px, y, py, \delta, \beta c\tau]$);
* ``` element_pass() ```Realiza o tracking de uma partícula um elemento do accelerador (objeto ```Element```);
*  ``` line_pass() ``` e ``` ring_pass() ``` Realiza o tracking de uma partícula em um acelerador (objeto ``` Accelerator ```);
* ``` find_orbit4() ``` e ``` find_orbit6() ``` Calcula a órbita fechada 4-D ou 6-D dos elétrons no acelerador.

> Logo mais serão adicionadas as funções de Tracking de multi-partículas em GPU (CUDA).

> Ainda estou escrevendo a documentação completa do pacote.
