#!/usr/local/bin/python

##
# Project:      pHOG
# Module:       DatabaseOp
# Desc.:        This module imports data stored in the db, creates tables, inserts
#                       entries according to the passed values
# Date: 18.02,02
# Modification history
# 19.02,02
# There are several issues not addressed:
#       1. No real error handeling mechanisms
#       2. The table creation method has no idea if the table exisited
#       3. The insert method doesn't check if the record intended for insertion is
#               redundant.
# 10.03.02
# Implement the second connect method, reading configuration from a file.
# 25.03.02
# Deal with password input from config file
##

import pg;
import _pg;
import sys;

# GLOBAL
connection = "";


##
# Read configuration file and return a dictionary. The ending "\n" is omitted.
##
def get_config(configFile,operation):

        configLines = open(configFile,"r").readlines()

        config = {}
        inOperation = 0
        for i in configLines:
                #allow documentation and empty string, match operation
                #print operation,i[:-1]
                if i[0] == ">":
                        if i[1:-1] != operation:
                                if inOperation:
                                        break
                                else:
                                        continue
                        elif i[1:-1] == operation:
                                inOperation = 1
                                continue

                #assign values to dict
                if inOperation and i[0] != "#" and i != "" and i != "\n":
                        if i[-1] == "\n":
                                config[i[:i.find(" ")]] = i[i.find(" ")+1:-1]
                        else:
                                config[i[:i.find(" ")]] = i[i.find(" ")+1:]

        return config



##
# Establish connection, need 4 parameters for it to work.
#
# @param  config - the dictionary with connection configuration
##
def connect(config):
        global connection;
        print "Establish connection..."
        
        connection = pg.connect(dbname = config["dbname"],
                                                        host   = config["host"],
                                                        port   = config["port"],
                                                        user   = config["user"])
        print "   Host = "+ connection.host
        print "     DB = "+ connection.db
        #print "   Port = %i" % connection.port
        #print " Status = %i\n" % connection.status
        return connection;

##
# In some situations, only dbname is passed.
##
def qconnect(config):

        global connection
        print "Establish quick connection..."

        connection = pg.connect(dbname = config["dbname"])
        print " Host    :",connection.host    
        print " Database:",connection.db

        return connection


##
# Establish connection. The config is in a file instead of a dictionary passed
# like above. The config file contains <config_name> <space> <value>
#
# @param configFile - the file with configuration information
# @param operation - read config from a particular operation
# @return config - the config dictionary
##
def configConnect(configFile, operation):
        global connection
        print "Establish connection.."
        configLines = open(configFile,"r").readlines()

        #store lines into a dict
        config = {}
        inOperation = 0
        for i in configLines:
                #allow documentation and empty string, match operation
                #print operation,i[:-1]
                if i[0] == ">":
                        if i[1:-1] != operation:
                                if inOperation:
                                        break
                                else:
                                        continue
                        elif i[1:-1] == operation:
                                inOperation = 1
                                continue

                #assign values to dict
                if inOperation and i[0] != "#" and i != "" and i != "\n":
                        if i[-1] == "\n":
                                config[i[:i.find(" ")]] = i[i.find(" ")+1:-1]
                        else:
                                config[i[:i.find(" ")]] = i[i.find(" ")+1:]

        if config.has_key("pass"):
                connection = pg.connect(dbname = config["dbname"],
                                                                host   = config["host"],
                                                                port   = int(config["port"]),
                                                                user   = config["user"],
                                                                passwd = config["pass"])
        else:
                connection = pg.connect(dbname = config["dbname"],
                                                                host   = config["host"],
                                                                port   = int(config["port"]),
                                                                user   = config["user"])
                
        print "   Host = "+ connection.host
        print "     DB = "+ connection.db
        #print "   Port = %i" % connection.port
        #print " Status = %i\n" % connection.status

        return config

##
# Close connection
##
def close():
        connection.close();


##
# Will send query to database according the passed argument. The sql command
# connect take the following arguments
# connect([dbname], [host], [port], [opt], [tty], [user], [passwd])
#
# @param  queryStr - the query string
# @return qTuple - the tuple returned from database query 
##
def select(queryStr,flag=0):

        if flag:
                print queryStr
        qTuple = (connection.query(queryStr)).getresult()

        return qTuple


def select2(fields,table,flag=0):

        queryStr = "SELECT %s FROM %s" % (fields,table)
        if flag:
                print queryStr
        qTuple = (connection.query(queryStr)).getresult()

        return qTuple


##
# Send query
##
def query(queryStr):
        #print "'%s'" % queryStr
        connection.query(queryStr)

##
# Creates a database with the specified name and fields
#
# @param tbName - name of the table to be created
# @param fields - a nested list contains the attributes and corresponding types
##
def createTable(tbName, fields):
        
        queryStr = "CREATE TABLE %s (" % tbName
        for i in fields:
                for j in i:
                        queryStr = queryStr+j+" "
                queryStr = queryStr[:-1]+","
        queryStr = queryStr[:-1]+");"
        print queryStr
        try:
                connection.query(queryStr)
                print "%s created successfully\n" % tbName
                return 1 #success
        except:
                print "%s already exist\n" % tbName
                
        

##
# Adding multiple entries to a database
# @param dbName name of the database to be inserted
#
# @param attrs  a string contains the attributes in the same order as the
#                               values to be inserted.
# @param values a nested list contains the attribute values, each element
#                               contains one entry (with multiple values) to be inserted
# @param flag   print debug string or not
##
def insert(tbName, attrs, values, flag = 0):
        for i in values:
                if flag:
                        print i
                queryStr = "INSERT INTO %s (%s) VALUES (" % (tbName, attrs)
                
                for j in i:
                        queryStr = queryStr+j+",";
                queryStr = queryStr[:-1]+");";

                if flag:
                        print " ",queryStr
                connection.query(queryStr);
                                
        return 1 #success


##
# Not tested yet
##
def dump(table):

        qTuple = self.dbtask.select("SELECT * FROM %s"% self.config[table])
        oup = open(self.config[table]+".dump","w")

        for i in qTuple:
                out_str = ""
                for j in i:
                        if type(j) == type(1):
                                j = str(j)
                        out_str = out_str + j + "\t"
                
                oup.write(out_str[:-1]+"\n")
                        
