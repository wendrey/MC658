// Resolvendo um programa linear pela biblioteca LEMON e GLPK
// Sugestões para Flávio K. Miyazawa, fkm@ic.unicamp.br  
// 
// Necessário instalar o LEMON e o GLPK (ou outro solver compatível). 
// O solver default é o GLPK. Depois,
// compile como: g++ ex_lp1_lemon.cpp -lemon -lglpk -o executavel
//
#include <lemon/lp.h>
using namespace lemon;
using namespace std;
int main()
{
  // Declarando/criando um LP (inicialmente vazio)
  Mip lp;
  string s;
  // Adicionando duas variáveis (colunas do sistema)
  Lp::Col x1 = lp.addCol();
  Lp::Col x2 = lp.addCol();
  lp.colType(x2, Mip::INTEGER);
  // Adicionando restrições ao programa linear
  lp.addRow(   x1 +    x2 <= 5.9 );
  lp.addRow( 2*x1 + 3*x2  <= 12.9);
  // Definindo os limitantes superior e inferior de cada variável
  lp.colLowerBound( x1, 0 );   lp.colUpperBound( x1, 4 );
  lp.colLowerBound( x2, 0 );
  // Definindo a função objetivo e se de minimização ou maximização
  lp.obj( -2*x1 - x2 );   lp.min();
  // Resolvendo o sistema
  lp.solve();

  // Imprimindo o valor das variáveis
  if (lp.type() == Mip::OPTIMAL) {
    cout << "Valor da funcao objetivo: " << lp.solValue() << endl;
    cout << "x1 = " << lp.sol(x1) << endl;
    cout << "x2 = " << lp.sol(x2) << endl;
  } else {cout << "Nao encontrou solucao otima." << endl;}
  return 0;
}
