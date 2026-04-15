// C implementation of
// the Quadratic Probing
#include <stdio.h>


// Function to print an array
void printArray(int arr[], int n)
{
	int i;
	// Iterating and printing the array
	for (i = 0; i < n; i++)
	{
		printf("%d\n", arr[i]);
	}
}

// Function to implement the
// quadratic probing
void hashing(int table[], int tsize,
			int arr[], int N)
{
	int i, j, hv, t;
	// Iterating through the array
	for (i = 0; i < N; i++)
	{
		// Computing the hash value
		hv = arr[i] % tsize;

		// Insert in the table if there
		// is no collision
		if (table[hv] == -1)
			{
				table[hv] = arr[i];
				printf("%d inserted at %d!\n", arr[i], hv);
			}
			
		else
		{
			// If there is a collision
			// iterating through all
			// possible quadratic values
			printf("Collision for %d at %d!\n", arr[i], hv);
			
			for (j = 0; j < tsize; j++)
			{
				// Computing the new hash value
				t = (hv + j * j) % tsize;
				if (table[t] == -1)
				{
					// Break the loop after
					// inserting the value
					// in the table
					table[t] = arr[i];
					printf("%d inserted at %d!\n", arr[i], t);
					break;
				}
				
				printf("Collision for %d at %d!\n", arr[i], t);
			}
		}
	}
	printf("After all keys are inserted, the table:\n\n");
	printArray(table, tsize);
}


// Function to implement search using
// quadratic probing
void search(int table[], int tsize,
			int key, int N)
{
	int i, j, hv, t, found = 0;
	
	
	// Computing the hash value
	hv = key % tsize;

	printf("\n\nSearch for key %d\n", key);
	
	if (table[hv] == key)
		{
			
			printf("%d Found at %d!\n", key, hv);
			found = 1;
		
		}
			
		else //collision
		{
			
			printf("After collision at %d!\n", hv);
			 
			for (j = 0; j < tsize; j++)
			{
			        
				// Computing the new hash value
				t = (hv + j * j) % tsize;
				if (table[t] == key)
				{
					
					printf("%d found at %d!\n", key, t);
					found = 1;
					break;
				}
				if(j>0){
				  printf("After collision at %d!\n", t);
				  }
				
			
			}
		}
		
		if(found == 0)
		{
			printf("Search failed!\n\n");
		}

}


//Deletion using quadratic probing
void deletion(int table[], int tsize,
			int key, int N)
{
	int i, j, hv, t, found = 0;

	
	// Computing the hash value
	hv = key % tsize;

	
	if (table[hv] == key)
		{
			table[hv] = -2; //Flag for deletion
			printf("%d deleted at %d!\n", key, hv);
			found = 1;
		
		}
			
		else //collision
		{
			
			
			for (j = 0; j < tsize; j++)
			{
				// Computing the new hash value
				t = (hv + j * j) % tsize;
				if (table[t] == key)
				{
					table[t] = -2; //flag for deletion
					printf("%d deleted at %d!\n", key, t);
					found = 1;
					break;
				}
				
			
			}
		}
		
		if(found == 0)
		{
			printf("Search failed!\n\n");
		}

}


int main()
{
	int arr[] = {46, 28, 21, 35, 57, 39, 19, 50};
	int N = 8, i, L = 11;

	// Size of the hash table
	
	int hash_table[11];

	// Initializing the hash table
	for (i = 0; i < L; i++)
	{
		hash_table[i] = -1;
	}

	// Quadratic probing
	hashing(hash_table, L, arr, N);
	
	search(hash_table, L, 57, N);
	
	//deletion(hash_table, L, 50, N);
	
	//printf("After deletion:\n\n");
	printArray(hash_table, L);
	return 0;
}



