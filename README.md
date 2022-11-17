Djeser Kordon 
Wei Liu

# Alignment free - TP 1

L'objectif du TP est de comparer 5 especes de bactéries entre elles.
Vous trouverez les données en suivant [ce lien](https://we.tl/t-WeWvheBBGX)

## Composer le TP

Vous devez forker ce projet puis compléter ses fonctions.
Le rendu sera le dépot git dans lequel vous aurrez forké.

Le but est d'obtenir toutes les distances paire à paire des différentes bactéries.
Vous pouvez modifier l'affichage final pour obtenir une matrice d'adjacence si vous les souhaitez.

En observant les distances obtenues, que pouvez-vous dire des espèces présentes dans cet échantillon ?

## Remarque sur le TME
 
Nous avons un problème dans le TME. Le premier est que la fonction est que l'on utilise ne compte pas les doublons ce qui fait qu'une valeur qui va apparaître plusieurs fois dans les deux listes fois va être compté Une seule fois dans le inter. Nous avons essayé une autre méthode (v2) pour calculer le inter mais elles ont des temps d'exécutions plus longs. Sachant que les différences sont inférieures à 1%.
 
Nous avons essayé de voir si d'autre méthode faisait que l'on avait des résultats différents et nous sommes arrivés à la conclusion que c'est le calcul du inter qui crée les différences de résultats. La deuxième méthode n'utilisant pas la fonction native de python, nous aurions tendance à plus lui faire confiance mais l'execution est plus longue.

ps: En fonction de l'OS, les fichiers ne sortent pas dans le même ordre. 
 
## interpretation des resultats
 
nous ne savons pas trop quoi dire des résultats. On a aucune information sur les espèces ainsi que leurs liens de parentés. On peut émettre des hypothèses en se basant sur les résultats mais on pourra juste les citer.Le coefficient (ou indice) de similarité de Jaccard est utilisé pour comparer la similarité entre deux ensembles donc dans notre cas, les deux ensemble de Kmer. et le coefficient de similitude la quantité d'un des ensembles dans le inter
