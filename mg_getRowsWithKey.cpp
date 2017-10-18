#include "mex.h"
#include "matrix.h"
#include <map>
#include <queue>
#include <string>

//  Author: Matthew Gunn
//  Updated:    11/10/15
//
//  Usage:      [key_unique, keymap] = mg_getRowsWithKey(keyvals)
//
//             This function 
//
//  keyvals:   This can either be (i)  a vector of real doubles
//                             or (ii) a cell array of strings
//
//  key_unique: This contains only the UNIQUE members of keyvals
//
//  keymap:     A cell array where keymap{i} will hold a vector containing
//              all the rows of keyvals that are equivalent to key_unique(i)
//
//  Example usage:
//    [key_unique, keymap] = mg_getRowsWithKey(group_id_vector)
// 
//    for(i=1:length(key_unique)),
//       group_i_key  = key_unique(i);
//       group_i_data = data(keymap{i},:);
//    end

void doitDouble(mxArray *plhs[], const mxArray *input) {

  size_t input_rows = mxGetM(input);

  double *input_data   = (double *) mxGetData(input);
  std::map<double, std::queue<size_t> > inputmap;                //map for output
  for(size_t i=0; i < input_rows; i++) {                         //iterate over the data, and push the row nums on
    if(mxIsNaN(input_data[i])) {
      mexErrMsgTxt("NaN was detected in input data. NaN keyvals aren't supported (currently would produce problematic behavior)");
    }
      
    inputmap[input_data[i]].push(i+1);
  }
  
  //create keyvals output
  size_t n_keys = inputmap.size(); 

  mxArray *output_array = mxCreateCellMatrix(n_keys, 1);                              // final cell array output

  mxArray *key_unique_array = mxCreateDoubleMatrix(n_keys, 1, mxREAL);  
  double *key_unique_data = mxGetPr(key_unique_array);

  std::map<double, std::queue<size_t> >::iterator imap_it;

  size_t j = 0;
  for(imap_it = inputmap.begin(); imap_it != inputmap.end(); imap_it++) {
    key_unique_data[j] = imap_it->first;           //assign the key unique value
    
    // ******************  In this section load up the vector associated with jth keyval ********
    std::queue<size_t> q = imap_it->second;   
    size_t n = q.size();   

    mxArray *group_j_array = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *group_j_data = mxGetPr(group_j_array);

    //  mexPrintf("%d\n",n);
    size_t k = 0;
    while(!q.empty()) {
      group_j_data[k] = q.front();
      q.pop();
      k++;
    }
    mxSetCell(output_array, j, group_j_array);             //set the cell

    j++;
  }

  plhs[0] = key_unique_array;
  plhs[1] = output_array;
}

void doitCellString(mxArray *plhs[], const mxArray *input) {

  size_t input_rows = mxGetM(input);
  
  //  double *input_data   = (double *) mxGetData(input);

  std::queue<char *> free_me_at_end;

  std::map<std::string, std::queue<size_t> > inputmap;                   //map for output
  for(size_t i=0; i < input_rows; i++) {                         //iterate over the data, and push the row nums on
    mxArray *rowi = mxGetCell(input, i);
    if(!mxIsChar(rowi))
      mexErrMsgTxt("Should be a cell array of strings");
    
    char *tempstr = mxArrayToString(rowi);                      //Allocates memory, we must call mxFree later
    free_me_at_end.push(tempstr);                               //so we can free later
    
    std::string asdf(tempstr);

    inputmap[asdf].push(i+1);    //pop it in
  }
  
  //create keyvals output
  size_t n_keys = inputmap.size(); 
  mxArray *output_array     = mxCreateCellMatrix(n_keys, 1);            //final cell array output
  mxArray *key_unique_array = mxCreateCellMatrix(n_keys, 1); 
  
  //  double *key_unique = (double *) mxCalloc(n_keys, sizeof(double));

  std::map<std::string, std::queue<size_t> >::iterator imap_it;
  size_t j = 0;

  for(imap_it = inputmap.begin(); imap_it != inputmap.end(); imap_it++) {
    //    key_unique[j] = imap_it->first;           //assign the key unique value

    mxArray *key_string_array = mxCreateString(imap_it->first.c_str());    

    mxSetCell(key_unique_array, j, key_string_array);   // FIXME: should I free key_string after this? (pretty sure no)

    // ******************  In this section load up the vector associated with jth keyval ********
    std::queue<size_t> q = imap_it->second;   
    size_t n = q.size();   
    //  mexPrintf("%d\n",n);

    mxArray *array   = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *row_data = mxGetPr(array);
    //    double *row_data = (double *) mxCalloc(n, sizeof(double));    

    size_t k = 0;
    while(!q.empty()) {
      row_data[k] = q.front();
      q.pop();
      k++;
    }
    //    mxSetData(array, row_data); 
    // *******************************************************************************************

    mxSetCell(output_array, j, array);             //FIXME: should I free after this? (i don't think so)
    j++;
  }
  plhs[0] = key_unique_array;
  plhs[1] = output_array;

  //free up memory
  while(!free_me_at_end.empty()) {
    mxFree(free_me_at_end.front());
    free_me_at_end.pop();
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if(nrhs != 1)
    mexErrMsgTxt("Invalid number of inputs.  Shoudl be 1 input argument.");
  
  if(nlhs != 2)
    mexErrMsgTxt("Invalid number of outputs.  Should be 2 output arguments.");
 
  const mxArray *input    = prhs[0];

  size_t input_cols = mxGetN(input);

  if(input_cols != 1)
    mexErrMsgTxt("Input should be 1 column");
  if(mxIsDouble(input)) {
    doitDouble(plhs, input);
  } else if(mxIsCell(input)) {
    doitCellString(plhs, input);
  } else {
    mexErrMsgTxt("Input should be double vector or cell array of strings");
  }
}











  
  /*
  double *output_data = (double *) mxCalloc(input_rows, sizeof(double));



  
  plhs[0] = 



  std::map<int, int>::iterator myit;

  int result;
  for(int i=0; i < input_rows; i++){
    myit = mymap.find(input_data[i]);
    if(myit == mymap.end())
      output_data[i] = 0;
    else
      output_data[i] = myit->second;
  }
  
  plhs[0] = mxCreateDoubleMatrix(input_rows, 1, mxREAL);
  mxSetData(plhs[0], output_data);*/
