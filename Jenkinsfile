pipeline {
    agent any

    stages {
        stage('Checkout') { // Checkout (git clone ...) the projects repository
            steps {
            checkout scm
            }
        }
        stage('Build') {
            steps {
                echo 'Building..'
                sh '''
                python -m venv env
                source env/bin/activate
                python -m pip install -r requirements.txt
                python -m pip install pylint
                '''
            }
        }

        stage('Test: Run') {
            steps {
                // Run Pylint.
                sh '''
                source env/bin/activate
                python -m pylint --exit-zero --format=parseable -d W1202 --ignore=config,deprecated,docs,env,tests,VaSeEval,VaSeUtils,vaseutils.py ./ > reports/pylint.report
                '''
            }
            post {
                always{
                    // Generate JUnit, PEP8, Pylint and Coverage reports.
                    recordIssues(
                        tools: [pyLint(pattern: 'reports/pylint.report')],
                        unstableTotalAll: 20,
                        failedTotalAll: 30
                    )
                }
            }
        }
        stage('Pylint') {
            steps {
                echo 'Testing..'
                sh '''
                source env/bin/activate
                python vase.py -h
                '''
                // python -m pylint --exit-zero -d W1202 --ignore=config,deprecated,docs,env,tests,VaSeEval,VaSeUtils,vaseutils.py
                // python -m pylint --ignore=config,deprecated,docs,env,tests,VaSeEval,VaSeUtils,vaseutils.py -d W1202,E1101 VaSeBuilder/
            }
        }
    }
}
