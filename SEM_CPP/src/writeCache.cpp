#include "iterativeSEM.hpp"
#include "../lib/sqlite3/sqlite3.h"
#include "common.hpp"
//#include <sstream>
using namespace std;

bool fileExists(const string &filename);
static void problemEncountered(const int &message, const string &what);
//static void isRowReady(const int &message);
// static void prepareStmt(sqlite3 *db, string stmt, sqlite3_stmt *query);
// static void checkDone(const int &message, const string &s);

// REQUIRES: accumSummary_scale is filled with the correct data
// EFFECTS: writes output of accumSummary_scale to cache, based upon dest
void writeCache(Dataset &data, const string &cache,
                Dataset::accumSummary_type::accumSummary_dest dest){

    // if(data.settings.verbose){
    //     cout << "Building cache for processed kmers.\n";
    // }
    // wants the third space, indexed from 0

    // points to an output vector from running accumSummary_scale(args)
    const vector<string> *ptr = nullptr;
    // const vector<double> *max_ptr = nullptr;
    switch (dest) {
        case Dataset::accumSummary_type::accumSummary_dest::none:
            cerr << "dest shouldn't be none!!!!" << endl;
            exit(1);
        break;
        case Dataset::accumSummary_type::accumSummary_dest::enumerated:
            ptr = &data.accumSummary_data.enum_accum_lines;
            // max_ptr = &data.accumSummary_data.enum_accum_max;
        break;
        case Dataset::accumSummary_type::accumSummary_dest::scrambled:
            ptr = &data.accumSummary_data.scramble_accum_lines;
            // max_ptr = &data.accumSummary_data.scramble_accum_max;
        break;
        case Dataset::accumSummary_type::accumSummary_dest::alignment:
            ptr = &data.accumSummary_data.align_accum_lines;
            // max_ptr = &data.accumSummary_data.align_accum_max;
        break;
        default:
            cerr << "there is no default for dest's switch statement!!!" << endl;
            exit(1);
        break;
    }
    // if(ptr->size() != max_ptr->size()){
    //     cerr << "ptr size mismatch!\n\tEXITING" << endl;
    //     switch(dest){
    //         case Dataset::accumSummary_type::accumSummary_dest::none:
    //         cerr << "none" << endl;
    //         exit(1);
    //     break;
    //     case Dataset::accumSummary_type::accumSummary_dest::enumerated:
    //         cerr << "enumerated" << endl;
    //     break;
    //     case Dataset::accumSummary_type::accumSummary_dest::scrambled:
    //         cerr << "scrambled" << endl;
    //     break;
    //     case Dataset::accumSummary_type::accumSummary_dest::alignment:
    //         cerr << "align" << endl;
    //     break;
    //     default:
    //         cerr << "default" << endl;;
    //         exit(1);
    //     break;
    //     }
    //     exit(1);
    // }

    // points to an output vector from running accumSummary_scale(args)

    bool newcache = fileExists(cache);

    sqlite3 *cacheDB;
    int message = 0;
    string msg;
    message = sqlite3_open(cache.c_str(), &cacheDB);
    problemEncountered(message, "open");

    //utilize sqlite transactions to speed this all up
    sqlite3_exec(cacheDB, "BEGIN TRANSACTION", NULL, NULL, NULL);
    sqlite3_exec(cacheDB, "PRAGMA journal_mode = MEMORY", NULL, NULL, NULL);

    if(newcache){

    }
    else{
        //We probably should never get to this build step because an initial cache is created earlier.
	// Are there any instances when this gets hit?

        if(data.settings.verbose){
            cout << "No existing cache and in WRITE CACHE(ERR)!\n";
        }
    // BLOB? why is it a BLOB? I used "text" not "blob", will need to ask about this

        // sqlite3_stmt *build_statement = nullptr;

        char *z_err_msg = NULL;

        msg = "CREATE TABLE kmer_cache (kmer TEXT PRIMARY KEY NOT NULL, alignment BLOB)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);

        msg = "CREATE UNIQUE INDEX kmerIDX ON kmer_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);

        msg = "CREATE TABLE seen_cache (kmer TEXT PRIMARY KEY NOT NULL, iter INT NOT NULL)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);

        msg = "CREATE UNIQUE INDEX seenIDX ON seen_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);
    }
    sqlite3_stmt *staged_query = nullptr;
    msg = "INSERT INTO kmer_cache VALUES(?,?)";
    sqlite3_prepare_v2(cacheDB, msg.c_str(), -1,
                       &staged_query, NULL);
    problemEncountered(message, msg);

    // Remove later per note below
    msg = "SELECT count(*) FROM kmer_cache WHERE kmer=? AND alignment=?";
    sqlite3_stmt *amount_seen_query = NULL;
    message = sqlite3_prepare_v2(cacheDB, msg.c_str(),
                                 static_cast<int>(msg.size()),
                                 &amount_seen_query, NULL);
    problemEncountered(message, msg);

    string temp = "", val = "";
    #ifdef DEBUG
    ofstream debug("written.txt");
    #endif

    for(size_t idx = 0; idx < ptr->size(); ++idx){

        val = ptr->at(idx);

        grab_string_4_index(val, temp);

        #ifdef DEBUG
        debug << "\tkmer: " << "first: #" << temp << "#\tsecond:#" << val << '#' << endl;
        #endif

        //currently we have kmers that have been seen in the input to this function
        // so these need to be filtered out here.
        // This should be edited later to not ever re-query these
        message = sqlite3_bind_text(amount_seen_query, 1, temp.c_str(),
                  -1, SQLITE_TRANSIENT);
        problemEncountered(message, "bind_text for amout_seen_query");

        message = sqlite3_bind_text(amount_seen_query, 2, val.c_str(),
                  -1, SQLITE_TRANSIENT);
        problemEncountered(message, "bind_text for amout_seen_query");
        message = sqlite3_step(amount_seen_query);

        #ifdef DEBUG
        if(sqlite3_column_type(amount_seen_query, 0) != SQLITE_INTEGER){
            cerr << "incorrect column_type on amount_seen_query\n\tEXITING"
                 << endl;
            exit(1);
         }
         #endif

        message = sqlite3_column_int(amount_seen_query, 0);

        if(message > 0) {
            //something found
        #ifdef DEBUG
            cerr << "DUP: " << message << endl;
        #endif
        } else {
            message = sqlite3_bind_text(staged_query, 1, temp.c_str(), -1, SQLITE_TRANSIENT);
            problemEncountered(message, "bind text 1 for staged_query, writeCache");
            message = sqlite3_bind_text(staged_query, 2, val.c_str(), -1, SQLITE_TRANSIENT);
            problemEncountered(message, "bind text 2 for staged_query, writeCache");
            message = sqlite3_step(staged_query);
            if(message != SQLITE_DONE){
                cerr << "Build statement is not done!\n\tEXITING" << endl;
                exit(1);
            }
            sqlite3_reset(staged_query);
            sqlite3_clear_bindings(staged_query);
        }

        sqlite3_reset(amount_seen_query);
        sqlite3_clear_bindings(amount_seen_query);
    }

    sqlite3_exec(cacheDB, "COMMIT TRANSACTION", NULL, NULL, NULL);

    message = sqlite3_finalize(staged_query);
    problemEncountered(message, "finalize staged_query");

    message = sqlite3_finalize(amount_seen_query);
    problemEncountered(message, "finalize amount_seen_query");

    message = sqlite3_close(cacheDB);
    problemEncountered(message, "closing the connection");
}



static void problemEncountered(const int &message, const string &what){
    if(message != SQLITE_OK){
        cerr << "Problem encountered with " << what << "!\n\tEXITING\n";
        cerr << "code: " << message << endl;
        exit(1);
    }
}

