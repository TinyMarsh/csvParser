#include <string>
#include <vector>

void getMascot(std::vector<std::string> &mascotSequences, std::vector<std::string> &mascotModPos, std::vector<std::string> &mascotCharges);
void getXTandem(std::vector<std::string> &xTandemSequences, std::vector<std::string> &xTandemModPos, std::vector<std::string> &xTandemCharges);

void getMascotQuants(std::vector<double> &mascotQuants);
void getXTandemQuants(std::vector<double> &xTandemQuants);

void getMascotRTs(std::vector<double> &mascotRTs);
void getXTandemRTs(std::vector<double> &xTandemRTs);

template<class Iter, class T>
Iter vectorBinarySearch(Iter begin, Iter end, T &val);