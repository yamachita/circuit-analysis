# Análise de Circuitos no Domínimo do Tempo

Programa para análise de circuitos no doíminio no tempo, aceita elementos lineares e não lineares. O sistema é implementado utilizando análise nodal modificada e o método dos trapézios junto com o método de Newton-Raphson, com técnica de *"Gmin stepping"* para ajudar a convergência.

## Elementos

- Fontes de corrente e tensão independentes.
- Fontes controladas.
- Resistores, capacitores e indutores.
- Amplificadores operacionais ideais, de 4 terminais.
- Transformadores ideais.
- Resistores lineares por partes.
- Chaves resistivas.

## Funcionamento

O sistema recebe como entrada um arquivo *.net* contendo uma *netlist* descrevendo o circuito, a saida do programa é um arquivo *.tab*, com o mesmo nome do arquivo de entrada, contendo a análise no tempo de cada elemento do circuito. No mesmo arquivo de entrada deve ser fornecido as instruções de como tratar a *netlist* conforme o formato explicado a seguir.

### Formato do Arquivo de Entrada

| Elemento                                    | Representação                                                                                                                 |
| ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| **Resistor**                                | `R<nome> <nó1> <nó2> <Resistência>`                                                                                           |
| **Indutor**                                 | `L<nome> <nó1> <nó2> <Indutância>`                                                                                            |
| **Capacitor**                               | `C<nome> <nó1> <nó2> <Capacitância>`                                                                                          |
| **Fonte de tensão controlada a tensão**     | `E<nome> <nóV+> <nóV-> <nóv+> <nóv-> <Av>`                                                                                    |
| **Fonte de corrente controlada a corrente** | `F<nome> <nóI+> <nóI-> <nói+> <nói-> <Ai>`                                                                                    |
| **Fonte de corrente controlada a tensão**   | `G<nome> <nóI+> <nóI-> <nóv+> <nóv-> <Gm>`                                                                                    |
| **Fonte de tensão controlada a corrente**   | `H<nome> <nóV+> <nóV-> <nói+> <nói-> <Rm>`                                                                                    |
| **Fonte de corrente**                       | `I<nome> <nó+> <nó-> <parâmetros>`                                                                                            |
| **Fonte de tensão**                         | `V<nome> <nó+> <nó-> <parâmetros>`                                                                                            |
| **Amplificador operacional ideal**          | `O<nome> <nó saída+> <nó saída-> <nó entrada+> <nó entrada->`                                                                 |
| **Resistor linear por partes**              | `N<nome> <nó+> <nó-> <4 pontos vi ji >`                                                                                       |
| **Transformador ideal**                     | `K<nome> <nó a> <nó b> <nó c> <nó d> <n>`                                                                                     |
| **Chave**                                   | `$<nome> <nó a> <nó b> <nó controle c> <nó controle d> <gon> <goff> <vref>`                                                   |
| **Fonte contínua**                          | `DC <valor>`                                                                                                                  |
| **Fonte senoidal**                          | `SIN <nível contínuo> <amplitude> <frequência em Hz> <atraso> <amortecimento> <defasagem em graus> <número de ciclos>`        |
| **Fonte pulsada**                           | `PULSE <amplitude 1> <amplitude 2> <atraso> <tempo de subida> <tempo de descida> <tempo ligada> <período> <número de ciclos>` |
| **Comentário**                              | `*<comentário>`                                                                                                               |
| **Configuração**                            | `.TRAN <tempo final> <passo> TRAP <passos por ponto na tabela>`                                                               |