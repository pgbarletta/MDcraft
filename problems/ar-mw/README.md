The steps to run the task:

1. **Prepare the data**
   Run the preparation script:
   ```bash
   python prepare.py
   ```

2. **Select and copy prepared atoms file**
   Choose one of the generated files matching the pattern prep*.h5 from the atoms/ folder and copy it to atoms.h5. For example:
   ```bash
   cp atoms/prep003000.h5 atoms.h5
   ```

3. **Run the main script**
   Execute the main program:
   ```bash
   python run.py
   ```
