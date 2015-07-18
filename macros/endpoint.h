// calculate the endpoint of two sequential 2-body decays of m2
float endpoint(float m2, float m1, float m0)
{
  return m2*sqrt( (1 - (m1/m2)*(m1/m2)) * (1 - (m0/m1)*(m0/m1)) );
}
