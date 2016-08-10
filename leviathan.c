#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Node_ {
   struct Node_ *parent;
   struct Node_ *left;
   struct Node_ *right;
   int isleaf;
   int arity;
   char content[100];
};
typedef struct Node_ node;

/*int printnode (node* n, char* str){
  char auxl[1000];
  char auxr[1000];
  char auxc[1000];
  if (n->isleaf){
    strcpy(str,n->content);
    return 1;
  }else{
    if(n->arity==1){
      strcpy(auxc,n->content);
      strcat(auxc,"(");
      printnode(n->left,auxl);printf("%s\n",auxl);
      strcat(auxc,auxl);
      strcat(auxc,")");
    }else{
      printnode(n->left,auxl);
      printnode(n->right,auxr);
      strcat(auxc,"(");
      strcat(auxc,auxl);
      strcat(auxc,n->content);
      strcat(auxc,")");
    }

  }
  strcpy(str,auxc);
  return 1;
}
*/

int printnode (node* n, char* str){
if(n->arity >= 1)
printnode(n->left, str);
strcat(str,n->content);
if (n->arity == 2)
printnode(n->right, str);
return 1;
}

void initnode(node* n,int height){
  //n= (node*) malloc(sizeof(node));
  strcpy(n->content,"a");
  if(height>0){
    n->left = (node*) malloc(sizeof(node));
    n->right = (node*) malloc(sizeof(node));
    (n->left)->parent = n;
    (n->right)->parent = n;
    initnode(n->left,height-1);
    initnode(n->right,height-1);
  }
}

int main() {
  char funcao[1000];
  strcpy(funcao,"");
  node *a;
  a = (node*) malloc(sizeof(node));
  initnode(a,4);
  strcpy(a->content,"*");
  a->arity =2;
  a->isleaf=0;
  strcpy((a->left)->content,"y");
  (a->left)->isleaf = 1;
  (a->left)->arity = 0;
  strcpy((a->right)->content,"exp");
  (a->right)->isleaf = 0;
  (a->right)->arity = 1;
  strcpy(((a->right)->left)->content,"*");
  ((a->right)->left)->isleaf=0;
  ((a->right)->left)->arity=2;
  strcpy((((a->right)->left)->left)->content,"-y");
  (((a->right)->left)->left)->isleaf = 1;
  (((a->right)->left)->left)->arity = 0;
  strcpy((((a->right)->left)->right)->content,"x");
  (((a->right)->left)->right)->isleaf=1;
  (((a->right)->left)->right)->arity=0;
  printnode(a,funcao);
  printf("%s\n",funcao);
  return 0;
//Hu3 BRBR
}

