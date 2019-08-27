main()
{
  int bit=1024;
  char *x;
  int first=1;
  double sum;
  
  while ((x)||(first==1)) 
  {
    first=0;
    x = malloc(bit);
    if (x)
    sum += bit;
  }
  printf("%20.10f bytes (%.1fMb)\n", sum, sum/1024.0/1024.0);
  return 0;
}

