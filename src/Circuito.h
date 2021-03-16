#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <string>
#include <vector>
#include <memory>
#include "Elemento.h"
using namespace std;

class Circuit
{

public:
    Circuit();
    void runSistema(string nome_arquivo);
    void printSolution();

    // Gets
    int getNumVariaveis() { return num_variaveis_; }
    int getNumElementos() { return num_elementos_; }
    int getNumNos() { return num_nos_; }

private:
    // Aquisicao de dados.
    vector<string> getFile(string nome_arquivo); // Pega a netlist de entrada do usuario.

    // Montagem NetList Interno
    void criarNetList(vector<string> arquivo); // Gera a netlist interna do sistema.
    void addCorrente();                        // Adiciona correntes como variaveis se necessario.
    void addElemento(char, const vector<string> &);
    int num(string); // Atribui numeros aos nós do sistema.

    // Montagem sistema
    void criarSistema(); // Cria o sistema matricial.
    void modSistemaTempo();
    void modSistemaNR();

    // Resolucao sistema
    bool resolverSistema(); // Resolve o sistema.

    void salvarSolucaoTempo(); // Salva a solucao atual na tabela de saida.
    void salvarTabela(int);
    void salvarSolucaoRN();

    // Auxiliar
    vector<string> copy(const vector<string> &, int, int);
    vector<string> split(string); // Separa uma string por espacos.
    string toUpperCase(string);

    // DEBUG
    void printNetList();
    void printSolucaoAtual(int i); // Printa o sistema resolvido no passo atual.
    void printSistemaInterno();    // Printa o NetList interno e as variaveis.

    // Modos
    void runDCLinear();
    void runDCNLinear();
    void runSPLinear();
    void runSPNLinear();
    double getErro(bool);

    // 1 - DC e Linear, 2 - DC e Não Linear
    // 3 - SIN/PULSE e Linear, 4 - SIN/PULSE e Não Linear
    int modo_;
    bool sp, linear;

    int num_variaveis_, num_elementos_, num_nos_;
    double dt_;
    double tempo_total_;
    double tempo_; // Tamanho do passo, tempo total, tempo atual.
    int passos_pontos_tabela;
    string nome_arquivo_;

    double gmin_;

    vector<unique_ptr<Elemento>> netlist1_;    // Elementos que a estampa não muda
    vector<unique_ptr<Elemento>> netlist2_;    // Elementos que a estampa muda com o tempo
    vector<unique_ptr<Elemento>> netlist3_;    // Elementos que a estampa muda no metodo de NR
    vector<vector<double>> sistema_;           // Contem o sistema a ser resolvido.
    vector<vector<double>> tabela_;            // Tabela com os valores das solucoes.
    vector<string> variaveis_;                 // Contem as variaveis do sistema.
    vector<vector<double>> sistema_resolvido_; // Contem a solução do sistema atual.
    vector<bool> convergiu_;

    vector<double> sol_atual_tempo_; // Contem a solução no tempo atual
    vector<double> sol_anterior_NR_; // Contem a solução anterior para o método de N-R
    vector<double> sol_atual_NR_;    // Contenm a solução atual para o método de N-R
};

#endif