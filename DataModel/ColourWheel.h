#ifndef COLOURWHEEL_H
#define COLOURWHEEL_H

#include <algorithm>
#include <TColor.h>
class ColourWheel{
	private:
	EColor currentcolour;
	int colourindex;
	static const std::vector<EColor> colours;
	static const std::vector<std::string> colournames;
	public:
	ColourWheel() : colourindex(0), currentcolour(colours.at(0)){};
	ColourWheel(int colourin){
		bool foundcolour = SetCurrentColour(colourin);
		if(foundcolour) return;
		colourindex=0;
		currentcolour = colours.at(colourindex);
	};
	ColourWheel(EColor colourin){
		bool foundcolour = SetCurrentColour(colourin);
		if(foundcolour) return;
		colourindex=0;
		currentcolour = colours.at(colourindex);
	};
	ColourWheel(std::string colourin){
		bool foundcolour = SetCurrentColour(colourin);
		if(foundcolour) return;
		colourindex=0;
		currentcolour = colours.at(colourindex);
	};
	~ColourWheel(){};
	void Reset(){ colourindex=1; }
	EColor GetNextColour(){
		colourindex++;
		if(colourindex>=colours.size()) colourindex=0;
		currentcolour = colours.at(colourindex);
		return currentcolour;
	};
	EColor GetPreviousColour(){
		colourindex--;
		if(colourindex<0) colourindex=colours.size()-1;
		currentcolour = colours.at(colourindex);
		return currentcolour;
	};
	EColor GetCurrentColour(){ return currentcolour; }
	int GetCurrentIndex(){ return colourindex; }
	bool SetCurrentColour(int index){
		if(index>0 && (index+1)<colours.size()){
			colourindex=index;
			currentcolour=colours.at(index);
			return true;
		} else { return false; }
	};
	bool SetCurrentColour(EColor colourin){
		auto colourpos = std::find(colours.begin(), colours.end(), colourin);
		if(colourpos!=colours.end()){
			colourindex=std::distance(colours.begin(),colourpos);
			currentcolour=colourin;
			return true;
		} else { return false; }
	};
	bool SetCurrentColour(std::string colourstring){
		auto colourpos = std::find(colournames.begin(), colournames.end(), colourstring);
		if(colourpos!=colournames.end()){
			colourindex=std::distance(colournames.begin(),colourpos);
			currentcolour=colours.at(colourindex);
			return true;
		} else { return false; }
	};
};

#endif
