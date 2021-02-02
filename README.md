**1.** **OverView**

In this project, we focus on using data structures and graph search algorithms to build a mapping application. We implemented several functions in this project, including auto complete location, get position, calculate shortest path between two nodes and travelling trojan problem.

Firstly, if the user doesn’t remember the exact name of a specific location, he or she can use the autocomplete function to find potential name of locations that start with or includes the input   characters for users to choose from. 

Secondly, if the user knows a location name, he or she can use the GetPosition function to get the latitude and longitude of this location. 

Thirdly, if users know two location’s name, he or she can use CalculateShortestPath function to get the shortest path between two locations. 

Lastly, users can input the number of locations and the TravellingTrojan function can give you the shortest route that covers all the locations and goes back to the start point and also the distance of this route(aka travelling salesman problem). 

We use a map to store special id and Node class which contains location information

This project is finished by DaSun and XinhuiTang. DaSun is responsible for GetPosition and TravellingTrojans function. Xinhui Tang is responsible for Autocomplete and CalculateShortestPath function. We finish our unit test by ourselves and combine them into our student test file. Also, we finished helper functions such as GetLat and GetLon etc. together So we think that we work on this project equally.

 

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image001.jpg)



 

 

 

**2.** **Classes and Helper functions**

**2.1 Node Class**

We use a class called Node to store the location information of one point in the map.

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image003.jpg)

Also, we use a map to store all the id and Node pair in our map in order to use id to visit all the information in the Node.

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image005.jpg)

 

 

 

**2.2 Helper functions:**

**1.GetLat:** In the data map, use location id to visit corresponding node’s latitude. 

**(Time complexity: O(1))**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image007.jpg)

 

**2. GetLon:** In the data map, use location id to visit corresponding node’s longitude. 

**(Time complexity: O(1))**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image009.jpg)

 

**3.GetName:** In the data map, use location id to visit corresponding node’s location name.

**(Time complexity: O(1))**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image011.jpg)

 

**4.GetNeighborIDs:** In the data map, use location id to acquire corresponding node’s neighbors. **(Time complexity: O(1))**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image013.jpg)

 

**5.CalculateDistance:** Use two node’s latitude and longitude and haversine formula to calculate the distance between two nodes. **(Time complexity: O(1))**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image015.jpg)

 

**6.CalculatePathLength:** traverse the input path vector and calculate the path length by adding the distance between every two nodes in the path. **(Time complexity: O(n))**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image017.jpg)

**3.** **Autocomplete**

**3.1 Autocomplete** **(****places that start with the search string)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image019.jpg)

In this Autocomplete function, we take the partial name of the location and return a list of possible locations with partial name as prefix. The principle is iterate data and find every location node, comparing it’s name’s prefix with partial name, if matches, add to results.

**Time complexity: O(n) (n is size of data)**

 

 

**Test Results (Autocomplete****(****places that start with the search string)):**

**1.ta**

![1606540747(1)](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image021.jpg)

 

**2.CA**

![1606540802(1)](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image023.jpg)

**3.2 Autocomplete2** **(****places that contain the search string)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image025.jpg)

I also extend the behavior so that it also shows the places that contains the search string, not just the ones start with the search string.

**Time complexity: O(n \* m) (n is the size of data, m is average length of location name)**

 

 

 

**Test Results (Autocomplete2** **(****places that contain the search string)):**

**1.fil**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image027.jpg)

 

**2.pa**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image029.jpg)

**4.** **GetPosition function** 

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image031.jpg)

In this GetPosition function, I use a pair of double to store the latitude and longitude. I use an iterator to traverse the data map to find the target location.

If the location is found, return the id corresponding to that location name and use function of GetLat and GetLon implemented before to set the output pair and return the pair.

If the location is not found, set the pair to be -1 and return the pair.

​    **Time complexity: O(n)**

 

 

 

 

 

 

 

 

**Test Results (GetPosition):**

**1.Target**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image033.jpg)

Target position in our map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image035.jpg)

Target position in google map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image037.jpg)

 

 

 

**2.ChickfilA**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image039.jpg)

ChickfilA position in our map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image041.jpg)

ChickfilA position in google map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image043.jpg)

 

 

 

**3.** **Popeyes Louisiana Kitchen**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image045.jpg)

Popeyes Louisiana Kitchen position in our map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image047.jpg)

Popeyes Louisiana Kitchen position in google map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image049.jpg)

 

 

**4. UCLA**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image051.jpg)

**5. CalculateShortestPath**

**5.1 CalculateShortestPath(Dijkstra):** 

/**

 \* CalculateShortestPath: Given 2 locations, return the shortest path which is a

 \* list of id. Dijkstra

 *

 \* @param {std::string} location1_name   : start

 \* @param {std::string} location2_name   : goal

 \* @return {std::vector<std::string>}    : path

 */

std::vector<std::string> TrojanMap::CalculateShortestPath(std::string location1_name, std::string location2_name) {

 std::vector<std::string> x;

 std::map<std::string, int> visited;

 std::map<std::string, double> dist;

 std::map<std::string, std::string> path;

 std::string srcID = GetID(location1_name);

 std::string destID = GetID(location2_name); 

 

 if(data.find(srcID) == data.end() || data.find(destID) == data.end() || srcID == "notFound" || destID == "notFound"){

  return x;

 }

 

 if(GetNeighborIDs(destID).size() == 0){

  return x;

 }

 

 for(auto it = data.begin(); it != data.end(); it++){

  dist[it->first] = DBL_MAX;

  visited[it->first] = 0;

 }

 

 visited[srcID] = 1;

 dist[srcID] = 0;

 

 std::vector<std::string> myNeighbors = GetNeighborIDs(srcID);

 for(auto n : myNeighbors){

  dist[n] = CalculateDistance(data[srcID], data[n]);

  path[n] = srcID;

 }

 

 for(auto it = dist.begin(); it != dist.end(); it++){

  std::string min_cost_id = findMinID(dist,visited);

  if(min_cost_id == destID){

   break;

  }

 

  std::vector<std::string> neighbors2 = GetNeighborIDs(min_cost_id);

  for(auto n : neighbors2){ 

   double distance = CalculateDistance(data[min_cost_id],data[n]);

   if(visited[n] != 1 && dist[min_cost_id] + distance < dist[n]){ 

​    dist[n] = dist[min_cost_id] + distance;

​    path[n] = min_cost_id;

   }

  } 

 }

 if(dist[destID] == DBL_MAX){

  return x;

 }

 x.push_back(destID);

 std::string t = path[destID];

 while(t != srcID){

  x.push_back(t);

  t = path[t];

 }

 x.push_back(srcID);

 

 std::reverse(x.begin(),x.end());

 

 return x;

}

 

// CalculateShortestPath helper function

std::string TrojanMap::GetID(std::string name){

 std::string ID;

std::string notFound = "notFound";

 if(name == ""){

  return notFound;

 }

 for(auto it = data.begin(); it != data.end(); it++){

  if(it->second.name == name){

   ID = it->first;

  }

 }

 return ID;

}

 

std::string TrojanMap::findMinID(std::map<std::string, double> &dist, std::map<std::string, int> &visited){

 double min = DBL_MAX;

 std::string minID;

 for(auto it = dist.begin(); it != dist.end(); it++){

 if(visited[it->first] != 1 && it->second < min){

 min = it->second;

 minID = it->first;

 }

 }

 visited[minID] = 1;

 return minID;

 

}

In this function, given 2 locations A and B, find the best route from A to B. My first solution uses Dijkstra algorithm. In this implementation, I set one map named visited to mark every extended node, one map named dist to record distance from source to every node in the map, one map named path to record node’s precursor node in the shortest path.

I also used a helper function called findMinID to find node which hasn’t been extended and has the shortest distance to source node, then mark node as extended. Return value is the ID of found node.

Dijkstra uses the strategy of greedy algorithm. In every time of iteration, one node which has the shortest distance(not extended yet) to source node is extended, then update it’s all neighbors nodes distance. When all nodes are extended, we can get shortest distance and path from source node to every node.

To decrease time consuming, once destination node is extended, algorithm can be stopped.

By find the precursor node of destination node recursively until we meet the start location, we can find the shortest path from start location to destination.

**Time Complexity: O(n^2)**

 

 

 

 

 

 

 

 

 

 

**Test Results (CalculateShortestPath(Dijkstra)):**

**1.Target--Ralphs**

The shortest path from Target to Ralphs in our map:

![1606541828(1)](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image053.jpg)![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image055.jpg)

The shortest path from Target to Ralphs in google map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image057.jpg)

 

 

 

 

**2.Cava—Target**

There is no path from Cava to Target using Dijkstra algorithm.

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image059.jpg)

 

 

 

 

 

 

 

 

 

 

 

 

 

 

**5. 2 CalculateShortestPath(Bellman-Ford):**

 

/**

 \* CalculateShortestPath: Given 2 locations, return the shortest path which is a

 \* list of id. Bellman_ford

 *

 \* @param {std::string} location1_name   : start

 \* @param {std::string} location2_name   : goal

 \* @return {std::vector<std::string>}    : path

 */

std::vector<std::string> TrojanMap::CalculateShortestPath_BF(

 std::string location1_name, std::string location2_name) {

 std::vector<std::string> x;

 std::map<std::string, double> dist;

 std::map<std::string, std::string> path;

 std::string srcID = GetID(location1_name);

 std::string destID = GetID(location2_name); 

if(data.find(srcID) == data.end() || data.find(destID) == data.end() || srcID == "notFound" || destID == "notFound"

){

  return x;

 }

 for(auto it = data.begin(); it != data.end(); it++){

  dist[it->first] = DBL_MAX;

 }

 dist[srcID] = 0;

 for(int i = 1; i <= dist.size()-1; i++ ){

  for(auto it = dist.begin(); it != dist.end(); it++){

   std::vector<std::string> myNeighbors = GetNeighborIDs(it->first);

   for(auto n : myNeighbors){

​    if( dist[n] > it->second + CalculateDistance(data[it->first], data[n])){

​     dist[n] = it->second + CalculateDistance(data[it->first], data[n]);

​     path[n] = it->first;

​    }

   }

  }

 }

 if(dist[destID] == DBL_MAX){

  return x;

 }

 x.push_back(destID);

 std::string t = path[destID];

 while(t != srcID){

  x.push_back(t);

  t = path[t];  }

 x.push_back(srcID);

 std::reverse(x.begin(),x.end());

 return x;

}

This is the second version of CalculateShortestPath using Bellman-ford algorithm. In this implementation, I use one map named dist to store distance from source node to other nodes and one map path to record precursor node. Dist should be initialized as dist[source] = 0 and dist[otherNode] = infinity.

Like Dijkstra’s algorithm, Bellman-Ford proceeds by relaxation, in which approximations to the correct distance are replaced by betters ones until they eventually reach the solution. This algorithm is much slower than Dijkstra, but it can handle situation that has negative edges. In this projects, every edge has positive weight so we can’t take this advantage. In my test, it takes about 40 seconds to reach final result while Dijkstra takes only about 5 seconds.

**Time Complexity:O(M\*N), M is number of vertex, N is number of edges**

 

 

 

 

 

 

 

 

 

 

 

 

**Test Results (CalculateShortestPath(Bellman-Ford)):**

**1.Popeyes Louisiana Kitchen--Target**

The shortest path from Popeyes Louisiana Kitchen to Target in our map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image061.jpg)![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image063.jpg)

 

 

The shortest path from Popeyes Louisiana Kitchen to Target in google map:

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image065.jpg)

 

 

 

**2.Cava—Ralphs**

There is no path from Cava to Target using Bellman-Ford algorithm.

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image067.jpg)

 

 

 

**5. 3 Comparison between Dijkstra and Bellman-Ford**

Comparison between Dijkstra and Bellman-Ford : Dijkstra algorithm’s nature is greedy algorithm, every time it extends node which has the shortest distance from source.So it can’t find shortest distance when there are negative edge in the graph. However, Bellman-Ford can handle this situation, it loops N-1 times(number of vertex) to do relaxation to each edge, if there is a negative cycle, after N-1 times loop, if distance can still be updated, it means there are negative cycle in the graph.However, in this project we didn’t use this feature.

Time complexity of Dijkstra is O(n^2), time complexity of Bellman-Ford is O(M*N), M is number of vertex, N is number of edges. Bellman-Ford consumes more time than Dijkstra. In this project, it takes about 4 seconds to find shortest path, but Bellman-Ford it takes about 45 seconds.

**6. TravellingTrojan**

**6. 1 TravellingTrojan (Brute Force with a helper function called TThelper)**

/**

 \* Travelling salesman problem: Given a list of locations, return the shortest

 \* path which visit all the places and back to the start point.

 *

 \* @param {std::vector<std::string>} input : a list of locations needs to visit

 \* @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path

 */

std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan(

​                  std::vector<std::string> &location_ids) {

 // Variables for TThelper

 std::pair<double, std::vector<std::vector<std::string>>> results;

 std::vector<std::vector<int>> resultsIndexVec;

 std::vector<int> path;

 std::vector<int> location_indexs; // list of locations of int

 std::vector<std::vector<double>> adjMatrix;

 

 // Create map between (string)id and (int)index to have access to adjacent matrix

 int map_index = 0;

 std::unordered_map<std::string, int> toNumMap; // string id to int

 for (auto id : location_ids){

  toNumMap[id] = map_index++;

 }

 

 // Using map to change input location_id to index

 for (auto id: location_ids){

  location_indexs.push_back(toNumMap[id]);

 }

 

 // Create map between (int)index and (string)id to have translate result path vector<int> into vector<string>

 int map_index2 = 0;

 std::unordered_map<int, std::string> toStrMap; // string id to int

 for (auto id : location_ids){

  toStrMap[map_index2++] = id;

 }

 

 // Create adjacent matrix according to location_ids

 // And using CalculateDistance function implemented before

 int size = location_ids.size();

 adjMatrix = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0)); // initialize adj martix to 0.0

 for(int i = 0; i < size; i++){

  for(int j = 0; j < size; j++){

   if(i != j){

​    // Create the distance value in adjacent matrix

​    adjMatrix[i][j] = CalculateDistance(data[location_ids[i]], data[location_ids[j]]); 

   } 

  }

 }

 

 // Using TSP helper to calculate 

 double minCost = DBL_MAX;

 double resLen = TTHelper(resultsIndexVec, location_indexs, adjMatrix, location_indexs[0], location_indexs[0], 0.0, minCost,path);

 

 // Transform int path into string path

 std::vector<std::vector<std::string>> resultsStr;

 std::vector<std::string> tmp;

 std::map<double, std::vector<int>>::iterator iter;

 

 // Transform int path into string path

 for(auto out: resultsIndexVec){

  for(auto inner: out){

   // Transform int path into string path

   tmp.push_back(toStrMap[inner]);

  }

  // Add one path into results

  resultsStr.push_back(tmp); 

  tmp.clear();

 }

 

 // Output results

 results = std::make_pair(resLen, resultsStr);

 

 return results;

} 

 

 

 

 

// TSP helper to calculate the Travelling Trojans Problem

double TrojanMap::TTHelper(std::vector<std::vector<int>> &resultsIndexVec, std::vector<int> &location_indexs, 

   std::vector<std::vector<double>> adjMatrix, int start, int curNode, double curCost, double &minCost, std::vector<int> path){

 

 double res = DBL_MAX;

 path.push_back(curNode);

 

 // If we have already visited all nodes, we add the last distance from last node to start node

 // then we compare current path length to minimum path length.

 // If curCost < minCost, we add start node the the last to current path to get a cuclic path

 //   and update the minimum path length

 if(path.size() == location_indexs.size()){

  double totalLen = curCost + adjMatrix[curNode][start];

  if(totalLen < minCost){

   std::vector<int> tmp(path);

   tmp.push_back(start);

   resultsIndexVec.push_back(tmp);

   minCost = totalLen;

  } 

  return totalLen; 

 }

 

 // iterate all nodes in current location index vector and generating branches

 for(int i = 0; i < location_indexs.size(); i++){

  // Get non-self and non-existing node in current path to construct rest path and using recursion to step into next node.

  // Avoid duplicate

  if(i != curNode && std::find(path.begin(), path.end(), i) == path.end()){

   // Pruning

   if(curCost < minCost){

​    // DFS like visiting all rest nodes and always add if next node is leaf node, return the shortest path to res variable. 

​    res = std::fmin(res, TTHelper(resultsIndexVec, location_indexs, adjMatrix, start, i, curCost + adjMatrix[curNode][i], minCost,path));

   }

  } 

 }

 

 return res;

}

 

The travelling salesman problem (also called the traveling salesperson problem or TSP) asks the following question: "Given a list of cities and the distances between each pair of cities, what is the shortest possible route that visits each city exactly once and returns to the origin city?" It is an NP-hard problem in combinatorial optimization.

In this function, I used Brute Force solution to solve the travelling salesman problem. Its key idea is evaluating all possible permutations and find the minimum.

I use an adjacent matrix called adjMatrix to store the distance data cause this data structure is very convenient when updating and getting the distance between every two nodes. Also, I used two maps called toNumMap and toStrMap to map (string)location id to (int)index and vice versa in order to use the adjMatrix convenitently.

In my helper function TThelper, I use **DFS** and **backtracking** to realize all possible paths and use **pruning** when finding that current cost have already bigger than minimum cost. When reaching the last node, accept the new path iff its cost is less than current cost. Return the shortest distance and along with all the paths when generating shorter paths.

**Time Complexity:O(n!)**

 

 

**Test** **Results(Brute Force):**

**1. Five Nodes(Brute Force)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image069.jpg)

The first picture shows the input number of places and the result of travelling salesman problem(5 nodes) in the form of location ids.

 

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image071.jpg)

The second picture is the result of the travelling salesman problem(5 nodes) in map view.

 

 

 

 

 

 

**2.Eight Nodes(Brute Force)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image073.jpg)

The first picture shows the input number of places and the result of travelling salesman problem(8 nodes) in the form of location ids.

 

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image075.jpg)

The second picture is the result of the travelling salesman problem(8 nodes) in map view.

 

 

 

 

 

**3.Ten Nodes(Brute Force)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image077.jpg)

​    The first picture shows the input number of places and the result of travelling salesman problem(10 nodes) in the form of location ids.

 

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image079.jpg)

​    The second picture is the result of the travelling salesman problem(10 nodes) in map view.

 

 

 

 

 

**6. 2 TravellingTrojan (2opt with a helper function called twoOptSwap)**

 

/**

 \* Travelling salesman problem: Given a list of locations, return the shortest

 \* path which visit all the places and back to the start point.

 *

 \* @param {std::vector<std::string>} input : a list of locations needs to visit

 \* @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path

 */

std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_2opt(

​                  std::vector<std::string> &location_ids) {

 std::pair<double, std::vector<std::vector<std::string>>> results;

 

 // Create a initial path which is the same order as input location_ids

 std::vector<std::string> curPath(location_ids);

 std::vector<std::string> newPath;

 std::vector<std::vector<std::string>> resultsPath;

 // put start node into path to build a cyclic route

 curPath.push_back(location_ids[0]);

 // Add initial path into result

 resultsPath.push_back(curPath);

 bool improveFlag = false;

 int size = curPath.size();

 double newCost;

 double bestCost = 0.0;

 // Using two pointer to iterate all pair of nodes except for the first and last one to do 2opt optimize.

 do{

  newCost = 0.0;

  improveFlag = false;

  bestCost = CalculatePathLength(curPath);

  for (int i = 1; i < size - 1; i++) {

​    for (int k = i + 1; k < size; k++) {

​      newPath = twoOptSwap(curPath, i, k);

​      newCost = CalculatePathLength(newPath);

​      // If we find a pair of switched node which leads to total distance shorter than current minimum distance,

​      // we add current path into result and update the current minimum distance.

​      // Also, we set two for loop to the end and update improveFlag to do another optimize loop

​      if (newCost < bestCost) {

​        curPath = newPath;

​        bestCost = newCost;

​        resultsPath.push_back(curPath);

​        improveFlag = true;

​        i = size - 1; // Goto do while again(Restart)

​        k = size;  // Goto do while again(Restart)

​      }

​    }

  }

 }while(improveFlag == true);

 

 // output results

 results = std::make_pair(bestCost, resultsPath);

 return results;

}

 

// swap path helper function

std::vector<std::string> TrojanMap::twoOptSwap(std::vector<std::string> curPath, int i, int k){

 std::vector<std::string> newPath(curPath);

 std::reverse(newPath.begin() + i, newPath.begin() + k);

 return newPath;

}

In this function, I used 2opt algorithm. The 2-opt algorithm basically removes two edges from the tour, and reconnects the two. You can see the picture of 2opt move below.

I use a vector to store my current path and use a vector<vector> to store the result. First, we generate an arbitrary path and iterate every two nodes and swap them to generate a new path. We compare the new path with current path. If the distance of new path is less than current path, update current path with new path and update current minimum cost with new cost, and then start from iterating every two nodes again until we have swapped all possible pair of nodes and there is no improve on the cost, we exit the algorithm.

This algorithm can just get **local minimum** instead of global optimum since it is **heuristic**

​    **Time Complexity:O(n^2)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image081.jpg)

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image083.jpg)![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image085.jpg)

**Test** **Results(2opt):**

**1.Five Nodes(2opt)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image087.jpg)

The first picture shows the input number of places and the result of travelling salesman problem(5 nodes) in the form of location ids.

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image089.jpg)

​    The second picture is the result of the travelling salesman problem(5 nodes) in map view.

 

 

 

 

 

 

**2.Ten Nodes(2opt)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image091.jpg)

The first picture shows the input number of places and the result of travelling salesman problem(10 nodes) in the form of location ids.

 

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image093.jpg)

​    The second picture is the result of the travelling salesman problem(10 nodes) in map view.

 

 

 

 

 

**3.Thirty Nodes(2opt)**

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image095.jpg)

The first picture shows the input number of places and the result of travelling salesman problem(30 nodes) in the form of location ids.

 

![img](file:////Users/sunda/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image097.jpg)

The second picture is the result of the travelling salesman problem(30 nodes) in map view.

**6. 3 Comparison between Brute Force and 2opt**

​    Time complexity: Brute Force takes O(n!) and 2opt takes O(n^2). There is a huge difference between these two algorithms because 2opt is just a heuristic and just swap different pair of nodes while brute force using DFS and backtracking to visit all possible paths. In practice, when we enter over 11 nodes, the brute force solution takes more than 10 seconds or even cannot solve the problem because of the huge time complexity while 2opt can handle even 30 to 40 nodes in less than 10 seconds.

​    As for the final result, 2opt can just give us the answer of local minimum and it doesn’t guarantee the result to be the global shortest path while brute force can give us global optimum, because brute force use DFS and backtracking and evaluate all possible permutations while 2opt is just heuristics which is not guaranteed to be optimal and is just used to generate some suboptimal solutions.

 

 

 

 

 

 

**7. Conclusion**

This project is finished by DaSun and XinhuiTang. DaSun is responsible for GetPosition and TravellingTrojans function. Xinhui Tang is responsible for Autocomplete and CalculateShortestPath function. We also finish our unit test by ourselves and combine them into our student test file. We think that we work on this project equally.

In this project, we learned several useful algorithms including Dijkstra, Bellman-Ford, Brute Force, 2opt, etc. We think that this project is a comprehensive practice for what we have learned this semester and is of great significance. When doing this project, we met some difficulties such as time cost to calculate shortest path is too long and sometimes it may ends with infinite loop. In order to solve those problems, we have to go through algorithm deeply and set break point to debug. Within the process of finding problem, solving problem, we deepen our understanding of shortest path algorithm and traveling salesman problem. Besides that, because work is divided, we need to write easy-understand code and annotation to help teammate understand. In this progress, we both learned how to work as a team efficiently.