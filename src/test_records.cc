#include <iostream>
#include <sys/stat.h>
#include "recordlocker.h"

#include "recordlocker.cc"

int main()
{
	mkdir("recordtest", 0777);

	recordlocker rl("recordtest");

	try
	{
		auto r = rl.lock_record();
		std::cout << "Got record " << r << std::endl;
		if(r)
		{
			sleep(20);
			rl.unlock_record();
		}
	}
	catch(std::string s)
	{
		std::cout << "error: " << s << "\n";
		return 1;
	}

	return 0;
}
