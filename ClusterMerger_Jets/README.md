# Jets Model
## Model 0
Consider constant current inside the injection region. Then the magnetic field is:
$$ \overrightarrow{B} = B_0 \frac{r}{R} \hat{\phi}, $$
where r is the radial distance, and R is the jets radius.\
Integrate through the injection region, the total magnetic energy has the form:
$$ E_B = \frac{B_0^2 R^3}{8}. $$

## Model 1 (not used)
Consider constant current inside A $(0 < A \leq 1)$ of the injection radius. (A is called `JetBModel1_AFactor` in the input file.) Then the magnetic field has the form:
$$ \overrightarrow{B} = \begin{cases} B_0 \frac{r}{a} \hat{\phi}, & \text{if r < a} \\ B_0 \frac{a}{r} \hat{\phi}, & \text{if r > a} \end{cases} $$
where $a = A * R$.\
Integrate through the injection region, the total magnetic energy has the form:
$$ E_B = \frac{R B_0^2}{2} (AR)^2 (\frac{1}{4} - ln(A)). $$