# Gatum

Gatum is a Java Based RESTful API for calculating average collision cross sections of molecules.

## Prerequisites

- Maven 3.5 or higher
- JDK 1.7

## Testing

This project could be run on "intellij idea community edition" by simply importing the "pom.xml" file to the IDE.
 
## Usage

A request to the API could be sent along with the chemical data file as follows:

curl -H "Content-Type: multipart/form-data" -F file=@"deprotonated_3D.sdf" http://localhost:8080/ccs

A sample data file could be found on the "input_file" folder
