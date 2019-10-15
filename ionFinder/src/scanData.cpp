//
//  scanData.cpp
//  ms2_anotator
//
//  Created by Aaron Maurais on 11/23/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#include <scanData.hpp>

void scanData::Scan::clear()
{
	_precursorFile.clear();
	_scanNum = 0;
	_sequence.clear();
	_charge = 0;
	_fullSequence.clear();
	_modified = false;
	_spectralCounts = 0;
}

std::string scanData::Scan::makeSequenceFromFullSequence(std::string fs) const
{
	fs = fs.substr(fs.find(".") + 1);
	fs = fs.substr(0, fs.find_last_of("."));
	return fs;
}

/**
 \brief Remove static modification symbols from peptide sequence \p s.
 
 \param s peptide sequence
 \param lowercase Should modified residue be transformed to lowercase?
 
 \return \s with static modifications removed.
 */
std::string scanData::removeStaticMod(std::string s, bool lowercase)
{
	//if s does not contain any modifications, just return s
	if(!(utils::strContains('(', s) ||
		 utils::strContains(')', s)))
		return s;
	
	std::string ret = "";
	for(int i = 0; i < s.length(); i++)
	{
		if(s[i] == ')')
			throw std::runtime_error("Invalid sequence: " + s);
		
		if(s[i] == '('){
			//get end of paren
			size_t end = s.find(')', i);
			if(end == std::string::npos)
				throw std::runtime_error("Invalid sequence: " + s);
			
			//erase paren from s
			s.erase(s.begin() + i, s.begin() + end + 1);
			
			if(lowercase)
				ret[ret.length() - 1] = std::tolower(ret.back());
		}
		ret += s[i];
	}
	
	return ret;
}

/**
 \brief Remove dynamic modification symbols from peptide sequence \p s.
 
 \param s peptide sequence
 \param lowercase Should modified residue be transformed to lowercase?
 
 \return \s with dynamic modifications removed.
 */
std::string scanData::removeDynamicMod(std::string s, bool lowercase)
{
	//if s does not contain any modifications, just return s
	if(!utils::strContains(scanData::MOD_CHAR, s))
		return s;
	
	std::string ret = "";
	for(int i = 0; i < s.length(); i++)
	{
		if(s[i] == scanData::MOD_CHAR){
			if(lowercase)
				ret[ret.length() - 1] = std::tolower(ret.back());
		}
		else ret += s[i];
	}
	
	return ret;
}

std::string scanData::Scan::makeOfSequenceFromSequence(std::string s) const{
	s = removeStaticMod(s);
	s = removeDynamicMod(s);
	return s;
}

void scanData::Scan::initilizeFromLine(std::string line)
{
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	_fullSequence = elems[12];
	_sequence = makeSequenceFromFullSequence(_fullSequence);
	_modified = utils::strContains(MOD_CHAR, _sequence);
	_xcorr = elems[2];
	_spectralCounts = std::stoi(elems[11]);
	
	std::string scanLine = elems[1];
	utils::split(scanLine, '.', elems);
	
	_precursorFile = elems[0] + ".ms2";
	_scanNum = std::stoi(elems[1]);
	_charge = std::stoi(elems[3]);
}

/**
 \brief Get unique output file name without extension.
 */
std::string scanData::Scan::getOfNameBase(std::string parentFile, std::string seq) const
{
	std::string ret;
	ret = utils::removeExtension(parentFile);
	ret += ("_" + makeOfSequenceFromSequence(seq) + "_" + std::to_string(_scanNum));
	if(_charge != 0)
		ret += ("_" + std::to_string(_charge));
	return ret;
}

/**
 \brief Calls getOfNameBase and adds OF_EXT to end.
 
 Only added for backwards compatibility.
 */
std::string scanData::Scan::getOfname() const{
	return getOfNameBase(_precursorFile, _sequence) + OF_EXT;
}
