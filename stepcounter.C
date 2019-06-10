void stepcounter(){
  double_t x = 0.018500;
  double_t mupdg = 4.*pow(0.1056583745,2);
  double_t inv_mass = sqrt(pow(x,2) + mupdg);

  while(inv_mass < 9.31){
    cout << x << " ";
    if(inv_mass < 1.0){inv_mass = inv_mass + 0.0005;}
    else {inv_mass = inv_mass + 0.001;}
  }
}
