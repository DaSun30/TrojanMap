#include "trojanmap.h"

#include <limits.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <locale>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <unordered_map>

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/videoio.hpp"

//-----------------------------------------------------
// TODO (Students): You do not and should not change the following functions:
//-----------------------------------------------------

/**
 * PrintMenu: Create the menu
 * 
 */
void TrojanMap::PrintMenu() {

  std::string menu =
      "**************************************************************\n"
      "* Select the function you want to execute.                    \n"
      "* 1. Autocomplete                                             \n"
      "* 2. Find the position                                        \n"
      "* 3. CalculateShortestPath                                    \n"
      "* 4. Travelling salesman problem                              \n"
      "* 5. Exit                                                     \n"
      "**************************************************************\n";
  std::cout << menu << std::endl;
  std::string input;
  getline(std::cin, input);
  char number = input[0];
  switch (number)
  {
  case '1':
  {
    menu =
        "**************************************************************\n"
        "* 1. Autocomplete                                             \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input a partial location:";
    std::cout << menu;
    getline(std::cin, input);
    auto results = Autocomplete(input);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.size() != 0) {
      for (auto x : results) std::cout << x << std::endl;
    } else {
      std::cout << "No matched locations." << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '2':
  {
    menu =
        "**************************************************************\n"
        "* 2. Find the position                                        \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input a location:";
    std::cout << menu;
    getline(std::cin, input);
    auto results = GetPosition(input);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.first != -1) {
      std::cout << "Latitude: " << results.first
                << " Longitude: " << results.second << std::endl;
      PlotPoint(results.first, results.second);
    } else {
      std::cout << "No matched locations." << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '3':
  {
    menu =
        "**************************************************************\n"
        "* 3. CalculateShortestPath                                            "
        "      \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input the start location:";
    std::cout << menu;
    std::string input1;
    getline(std::cin, input1);
    menu = "Please input the destination:";
    std::cout << menu;
    std::string input2;
    getline(std::cin, input2);
    auto results = CalculateShortestPath(input1, input2);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.size() != 0) {
      for (auto x : results) std::cout << x << std::endl;
      PlotPath(results);
    } else {
      std::cout << "No route from the start point to the destination."
                << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '4':
  {
    menu =
        "**************************************************************\n"
        "* 4. Travelling salesman problem                              \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "In this task, we will select N random points on the map and you need to find the path to travel these points and back to the start point.";
    std::cout << menu << std::endl << std::endl;
    menu = "Please input the number of the places:";
    std::cout << menu;
    getline(std::cin, input);
    int num = std::stoi(input);
    std::vector<std::string> keys;
    for (auto x : data) {
      keys.push_back(x.first);
    }
    std::vector<std::string> locations;
    srand(time(NULL));
    for (int i = 0; i < num; i++)
      locations.push_back(keys[rand() % keys.size()]);
    PlotPoints(locations);
    std::cout << "Calculating ..." << std::endl;
    auto results = TravellingTrojan(locations);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    CreateAnimation(results.second);
    if (results.second.size() != 0) {
      for (auto x : results.second[results.second.size()-1]) std::cout << x << std::endl;
      menu = "**************************************************************\n";
      std::cout << menu;
      std::cout << "The distance of the path is:" << results.first << std::endl;
      PlotPath(results.second[results.second.size()-1]);
    } else {
      std::cout << "The size of the path is 0" << std::endl;
    }
    menu = "**************************************************************\n"
           "You could find your animation at src/lib/output.avi.          \n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '5':
    break;
  default:
    std::cout << "Please select 1 - 5." << std::endl;
    PrintMenu();
    break;
  }
}


/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 * 
 */
void TrojanMap::CreateGraphFromCSVFile() {
  std::fstream fin;
  fin.open("src/lib/map.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '['), word.end());
      word.erase(std::remove(word.begin(), word.end(), ']'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}

/**
 * PlotPoint: Given a location id, plot the point on the map
 * 
 * @param  {std::string} id : location id
 */
void TrojanMap::PlotPoint(std::string id) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto result = GetPlotLocation(data[id].lat, data[id].lon);
  cv::circle(img, cv::Point(result.first, result.second), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}
/**
 * PlotPoint: Given a lat and a lon, plot the point on the map
 * 
 * @param  {double} lat : latitude
 * @param  {double} lon : longitude
 */
void TrojanMap::PlotPoint(double lat, double lon) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto result = GetPlotLocation(lat, lon);
  cv::circle(img, cv::Point(int(result.first), int(result.second)), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  cv::startWindowThread();
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}

/**
 * PlotPath: Given a vector of location ids draws the path (connects the points)
 * 
 * @param  {std::vector<std::string>} location_ids : path
 */
void TrojanMap::PlotPath(std::vector<std::string> &location_ids) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto start = GetPlotLocation(data[location_ids[0]].lat, data[location_ids[0]].lon);
  cv::circle(img, cv::Point(int(start.first), int(start.second)), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  for (auto i = 1; i < location_ids.size(); i++) {
    auto start = GetPlotLocation(data[location_ids[i - 1]].lat, data[location_ids[i - 1]].lon);
    auto end = GetPlotLocation(data[location_ids[i]].lat, data[location_ids[i]].lon);
    cv::circle(img, cv::Point(int(end.first), int(end.second)), DOT_SIZE,
               cv::Scalar(0, 0, 255), cv::FILLED);
    cv::line(img, cv::Point(int(start.first), int(start.second)),
             cv::Point(int(end.first), int(end.second)), cv::Scalar(0, 255, 0),
             LINE_WIDTH);
  }
  cv::startWindowThread();
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}

/**
 * PlotPoints: Given a vector of location ids draws the points on the map (no path).
 * 
 * @param  {std::vector<std::string>} location_ids : points
 */
void TrojanMap::PlotPoints(std::vector<std::string> &location_ids) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  for (auto x : location_ids) {
    auto result = GetPlotLocation(data[x].lat, data[x].lon);
    cv::circle(img, cv::Point(result.first, result.second), DOT_SIZE,
               cv::Scalar(0, 0, 255), cv::FILLED);
  }
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}


/**
 * CreateAnimation: Create the videos of the progress to get the path
 * 
 * @param  {std::vector<std::vector<std::string>>} path_progress : the progress to get the path
 */
void TrojanMap::CreateAnimation(std::vector<std::vector<std::string>> path_progress){
  cv::VideoWriter video("src/lib/output.avi", cv::VideoWriter::fourcc('M','J','P','G'), 10, cv::Size(1248,992));
  for(auto location_ids: path_progress) {
    std::string image_path = cv::samples::findFile("src/lib/input.jpg");
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    auto start = GetPlotLocation(data[location_ids[0]].lat, data[location_ids[0]].lon);
    cv::circle(img, cv::Point(int(start.first), int(start.second)), DOT_SIZE,
              cv::Scalar(0, 0, 255), cv::FILLED);
    for (auto i = 1; i < location_ids.size(); i++) {
      auto start = GetPlotLocation(data[location_ids[i - 1]].lat, data[location_ids[i - 1]].lon);
      auto end = GetPlotLocation(data[location_ids[i]].lat, data[location_ids[i]].lon);
      cv::circle(img, cv::Point(int(end.first), int(end.second)), DOT_SIZE,
                cv::Scalar(0, 0, 255), cv::FILLED);
      cv::line(img, cv::Point(int(start.first), int(start.second)),
              cv::Point(int(end.first), int(end.second)), cv::Scalar(0, 255, 0),
              LINE_WIDTH);
    }
    video.write(img);
    cv::imshow("TrojanMap", img);
    cv::waitKey(1);
  }
	video.release();
}
/**
 * GetPlotLocation: Transform the location to the position on the map
 * 
 * @param  {double} lat         : latitude 
 * @param  {double} lon         : longitude
 * @return {std::pair<double, double>}  : position on the map
 */
std::pair<double, double> TrojanMap::GetPlotLocation(double lat, double lon) {
  std::pair<double, double> bottomLeft(34.01001, -118.30000);
  std::pair<double, double> upperRight(34.03302, -118.26502);
  double h = upperRight.first - bottomLeft.first;
  double w = upperRight.second - bottomLeft.second;
  std::pair<double, double> result((lon - bottomLeft.second) / w * 1248,
                                   (1 - (lat - bottomLeft.first) / h) * 992);
  return result;
}

//-----------------------------------------------------
// TODO: Student should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id.
 * 
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(std::string id) {
  return data[id].lat;
}

/**
 * GetLon: Get the longitude of a Node given its id. 
 * 
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */
double TrojanMap::GetLon(std::string id) {
  return data[id].lon;
}

/**
 * GetName: Get the name of a Node given its id.
 * 
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(std::string id) {
  return data[id].name;
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node.
 * 
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(std::string id) {
    std::vector<std::string> result = data[id].neighbors;
    return result;
}


/**
 * CalculateDistance: Get the distance between 2 nodes. 
 * 
 * @param  {Node} a  : node a
 * @param  {Node} b  : node b
 * @return {double}  : distance in mile
 */
double TrojanMap::CalculateDistance(const Node &a, const Node &b) {
  // TODO: Use Haversine Formula:
  // dlon = lon2 - lon1;
  // dlat = lat2 - lat1;
  // a = (sin(dlat / 2)) ^ 2 + cos(lat1) * cos(lat2) * (sin(dlon / 2)) ^ 2;
  // c = 2 * arcsin(min(1, sqrt(a)));
  // distances = 3961 * c;

  // where 3961 is the approximate radius of the earth at the latitude of
  // Washington, D.C., in miles
  double dLon = (b.lon - a.lon)* M_PI / 180.0;
  double dLat = (b.lat - a.lat)* M_PI / 180.0;
  double lat1 = (a.lat) * M_PI / 180.0; 
  double lat2 = (b.lat) * M_PI / 180.0; 
  double d = std::pow(std::sin(dLat / 2),2) + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(dLon / 2),2);
  double c = 2 * std::asin(std::fmin(1, std::sqrt(d)));
  double distances = 3961 * c;
  return distances;
}

/**
 * CalculatePathLength: Calculates the total path length for the locations inside the vector.
 * 
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) {
  double len = 0.0;
  for(int i = 0; i < path.size() - 1; i++){
    len += CalculateDistance(data[path[i]], data[path[i+1]]);
  }
  return len;
}

/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix.
 *
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */
std::vector<std::string> TrojanMap::Autocomplete(std::string name) {
  std::vector<std::string> results;
  std::transform(name.begin(),name.end(),name.begin(), ::tolower);
  std::map<std::string, Node>::iterator it;
  it = data.begin();
  while(it != data.end()){
    Node location = it->second;
    std::string location_name = location.name;
    std::transform(location_name.begin(),location_name.end(),location_name.begin(), ::tolower);
    if(location_name.substr(0,name.length()) == name){
      results.push_back(location.name);
    }
    it++;
  }
  return results;
}


// AutoComplete2 can show the places that contains the search string, not just the ones that start with the search string
std::vector<std::string> TrojanMap::Autocomplete2(std::string name) {
  std::vector<std::string> results;
  int length = name.length();
  std::transform(name.begin(),name.end(),name.begin(), ::tolower);
  std::map<std::string, Node>::iterator it;
  std::map<std::string, Node>::iterator start;
  it = data.begin();
  while(it != data.end()){
    Node location = it->second;
    std::string location_name = location.name;
    std::transform(location_name.begin(),location_name.end(),location_name.begin(), ::tolower);
    if(location_name.length() >= length){
      for(int i = 0; i <= location_name.length() - length; i++){
        if(location_name.substr(i,name.length()) == name){
          results.push_back(location.name);
        }
      }
    }
    it++;
  }
  return results;
}




/**
 * GetPosition: Given a location name, return the position.
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::pair<double, double> results;
  std::string id_;
  std::map<std::string, Node>::iterator it;
  for(it = data.begin(); it != data.end(); it++){
    if(it->second.name == name){
      id_ = it->second.id;
      break;
    }
  }
  if (it == data.end()){
    results.first = -1;
    results.second = -1;
  }else{
    results.first = GetLat(id_);
    results.second = GetLon(id_);
  }
  return results;
}


/**
 * CalculateShortestPath: Given 2 locations, return the shortest path which is a
 * list of id.  Dijkstra
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath(std::string location1_name, std::string location2_name) {
  std::vector<std::string> x;
  std::map<std::string, int> visited;
  std::map<std::string, double> dist;
  std::map<std::string, std::string> path;
  std::string srcID = GetID(location1_name);
  std::string destID = GetID(location2_name); 
  
  if(data.find(srcID) == data.end() || data.find(destID) == data.end()){
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
        dist[n] = dist[min_cost_id] + distance;
        path[n] = min_cost_id;
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
  for(auto it = data.begin(); it != data.end(); it++){
    if(it->second.name == name){
      ID = it->first;
    }
  }
  return ID;
}

bool TrojanMap::isAllVisited(std::map<std::string, bool> visited){
  for(auto it = visited.begin(); it != visited.end(); it++){
    if(it->second == false)    return false;
  }
  return true;
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



/**
 * CalculateShortestPath: Given 2 locations, return the shortest path which is a
 * list of id. Bellman_ford
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_BF(
  std::string location1_name, std::string location2_name) {
  std::vector<std::string> x;
  std::map<std::string, double> dist;
  std::map<std::string, std::string> path;
  std::string srcID = GetID(location1_name);
  std::string destID = GetID(location2_name); 
  
  for(auto it = data.begin(); it != data.end(); it++){
    dist[it->first] = DBL_MAX;
  }
  dist[srcID] = 0;
  for(int i = 1; i <= dist.size()-1; i++ ){
  // std::map<std::string, std::string> visitedEdge = {};
  
    for(auto it = dist.begin(); it != dist.end(); it++){
      std::vector<std::string> myNeighbors = GetNeighborIDs(it->first);
      for(auto n : myNeighbors){
        if( dist[n] > it->second + CalculateDistance(data[it->first], data[n])){
          dist[n] = it->second + CalculateDistance(data[it->first], data[n]);
          path[n] = it->first;
        }
      }
    }
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





/**
 * Travelling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path
 */
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan(
                                    std::vector<std::string> &location_ids) {
  // Variables for TThelper
  std::pair<double, std::vector<std::vector<std::string>>> results;
  std::vector<std::vector<int>> resultsIndexVec;
  std::vector<int> path;
  std::vector<int> location_indexs;  // list of locations of int
  std::vector<std::vector<double>> adjMatrix;

  // Create map between (string)id and (int)index to have access to adjacent matrix
  int map_index = 0;
  std::unordered_map<std::string, int> toNumMap;  // string id to int
  for (auto id : location_ids){
    toNumMap[id] = map_index++;
  }

  // Using map to change input location_id to index
  for (auto id: location_ids){
    location_indexs.push_back(toNumMap[id]);
  }

  // Create map between (int)index and (string)id to have translate result path vector<int> into vector<string>
  int map_index2 = 0;
  std::unordered_map<int, std::string> toStrMap;  // string id to int
  for (auto id : location_ids){
    toStrMap[map_index2++] = id;
  }

  // Create adjacent matrix according to location_ids
  // And using CalculateDistance function implemented before
  int size = location_ids.size();
  adjMatrix = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));  // initialize adj martix to 0.0
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      if(i != j){
        // Create the distance value in adjacent matrix
        adjMatrix[i][j] = CalculateDistance(data[location_ids[i]], data[location_ids[j]]); 
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
  //     and update the minimum path length
  if(path.size() == location_indexs.size()){
    double totalLen = curCost + adjMatrix[curNode][start];
    if(totalLen < minCost){
      std::vector<int> tmp(path);
      tmp.push_back(start);
      resultsIndexVec.push_back(tmp);
      // 展示vec
      // std::cout << "个path为:[";
      // for(auto inner: tmp){
      //   // tmp.push_back(toStrMap[inner]); // 翻译一条路径
      //   std::cout << inner << " ->";
      // }
      minCost = totalLen;
    } 
    //std::cout << "resultIndexVec的size:" << resultsIndexVec.size() << std::endl;
    //std::cout << "totalLen:" << totalLen << std::endl;
    return totalLen; 
  }

  // iterate all nodes in current location index vector and generateing branches
  for(int i = 0; i < location_indexs.size(); i++){
    // Get non-self and non-existing node in current path to construct rest path and using recursion to step into next node.
    // Avoid duplicate
    if(i != curNode && std::find(path.begin(), path.end(), i) == path.end()){
      // Pruning
      if(curCost < minCost){
        // DFS like visiting all rest nodes and always add if next node is leaf node, return the shortest path to res variable.  
        res = std::fmin(res, TTHelper(resultsIndexVec, location_indexs, adjMatrix, start, i, curCost + adjMatrix[curNode][i], minCost,path));
      }
      
    } 
  }

  return res;
}




/**
 * Travelling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path
 */
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_2opt(
                                    std::vector<std::string> &location_ids) {
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
    // std::cout << "bestCost:" << bestCost << std::endl;
    // std::cout << "size:" << size << std::endl;
    for (int i = 1; i < size - 1; i++) {
        for (int k = i + 1; k < size; k++) {
            newPath = twoOptSwap(curPath, i, k);
            newCost = CalculatePathLength(newPath);
            // std::cout << "newCost:" << newCost;
            // If we find a pair of switched node which leads to total distance shorter than current minimum distance,
            // we add current path into result and update the current minimum distance.
            // Also, we set two for loop to the end and update improveFlag to do another optimize loop
            if (newCost < bestCost) {
                curPath = newPath;
                bestCost = newCost;
                resultsPath.push_back(curPath);
                improveFlag = true;
                i = size - 1;  // Goto do while again(Restart)
                k = size;   // Goto do while again(Restart)
            }
            // std::cout << "improveFlagInnerFor:" << improveFlag << std::endl;
        }
        // std::cout << "improveFlagOuterFor:" << improveFlag << std::endl;
    }
    // std::cout << "improveFlag:" << improveFlag << std::endl;
  }while(improveFlag == true);

  // output results
  results = std::make_pair(bestCost, resultsPath);

  return results;

}

// swap path helper function
std::vector<std::string> TrojanMap::twoOptSwap(std::vector<std::string> curPath, int i, int k){
  // std::cout << "    i:" << i << "k:" << k << std::endl;
  std::vector<std::string> newPath(curPath);
  std::reverse(newPath.begin() + i, newPath.begin() + k);
  //展示vec
  // std::cout << "个path为:[";
  // for(auto inner: newPath){
  //   std::cout << inner << " ->";
  // }
  return newPath;
}